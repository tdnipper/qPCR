library(tidyverse)

# Read in the data
data <- read_csv("TN045/45.8/TN045.8_filtered.csv")

# Summarize mean CT
data <- data|>
  group_by(`Sample Name`, `Target Name`, Task) |>
  summarize(mean_ct = mean(CT), .groups = "drop")

# Get non housekeeping targets
targets <- data |> 
  filter(`Target Name` != "TBP1") |>
  pull(`Target Name`) |>
  unique()

# Pivot data to have one column per target and calculate dCT by subtracting TBP1 from each target
data_dct <- data |>
  pivot_wider(names_from = `Target Name`, values_from = mean_ct) |>
  group_by(`Sample Name`, Task) |>
  mutate(
	across(
		.cols = all_of(targets), 
		.fn = ~ . - TBP1, 
		.names = "dCT_{col}")
		) |>
  ungroup()

# Calculate mock mean dCT of DUSP11 and PB2 and store as named vector
mock_values <- data_dct |>
  filter(`Sample Name` == "mock") |>
  summarize(
	across(
		.cols = starts_with("dCT_"), 
		mean)
		) |>
  unlist()

# Calculate ddCT by subtracting the dCT of the mock sample from each condition
data_ddct <- data_dct |>
  mutate(
	across(
		starts_with("dCT_"),
		~ . - mock_values[cur_column()],
		.names = "ddCT_{col}"
	)
  ) |>
  select(`Sample Name`, Task, starts_with("ddCT_")) |>
  rename_with(
	.cols = starts_with("ddCT_"),
	.fn = ~ str_replace(.x, "ddCT_dCT_", "ddCT_")
  )

# Calculate fold change as 2^(-ddCT)
data_foldchange <- data_ddct |>
  mutate(
	across(
		starts_with("ddCT_"),
		~ 2^(-.),
		.names = "foldchange_{col}"
	)
  ) |>
  select(`Sample Name`, Task, starts_with("foldchange_")) |>
  rename_with(
	.cols = starts_with("foldchange_ddCT_"),
	.fn = ~ str_replace(.x, "foldchange_ddCT_", "foldchange_")
  )

# Write ddCT and fold change data to csv
data_write <- data_ddct |>
  left_join(data_foldchange, by = c("Sample Name", "Task"))
write_csv(data_write, "TN045/45.8/TN045.8_ddCT_foldchange.csv")

# Plot fold change for DUSP11 and PB2	
p_DUSP11 <- ggplot(
	data_foldchange, 
	aes(
		x = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS", "del-NS1")),
		y = foldchange_DUSP11)) +
	geom_point(
		position = position_jitter(width = 0.1)) +
  labs(
		x = "", 
		y = "Fold Change") +
  theme_classic()


ggsave("TN045/45.8/TN045.8_DUSP11_foldchange_plot.png", plot = p_DUSP11, width = 8, height = 6, dpi = 300)


p_PB2 <- ggplot(
	data_foldchange, 
	aes(
		x = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS", "del-NS1")),
		y = foldchange_WSN_PB2)) +
	geom_point(
		position = position_jitter(width = 0.1)) +
	labs(
		x = "", 
		y = "Fold Change") +
	scale_y_log10() +
	theme_classic()

ggsave("TN045/45.8/TN045.8_PB2_foldchange_plot.png", plot = p_PB2, width = 8, height = 6, dpi = 300)

# Plot PB2 ddCT of PA-FS and del-NS1 relative to WT, not mock
# Get dCT values for new calculation
data_PB2_ratio <- data_dct |>
  filter(`Sample Name` != "mock") |>
	select(`Sample Name`, Task, dCT_WSN_PB2) |>
	pivot_wider(names_from = `Sample Name`, values_from = dCT_WSN_PB2)

# Calculate mean dCT of WT for PB2 to use as reference for ddCT calculation
PB2_WT_mean <- data_PB2_ratio |>
  pull(WT) |>
	mean()

# Calculate ddCT for PA-FS and del-NS1 relative to WT and then calculate fold change
data_PB2_ratio <- data_PB2_ratio |>
	mutate(
	PA_FS_ddCT = `PA-FS` - PB2_WT_mean,
	del_NS1_ddCT = `del-NS1` - PB2_WT_mean
	) |>
	select(Task, PA_FS_ddCT, del_NS1_ddCT) |>
	mutate(
		across(
			.cols = ends_with("ddCT"),
			.fn = ~ 2^(-.),
			.names = "foldchange_{col}"
		)) |>
	rename_with(
		.cols = starts_with("foldchange"),
		.fn = ~ str_replace(.x, "_ddCT", "")
	)

# Pivot longer for plotting
data_PB2_ratio_long <- data_PB2_ratio |>
	pivot_longer(
	cols = starts_with("foldchange"),
	names_to = "strain",
	values_to = "FoldChange"
	) |>
	mutate(
	strain = str_replace(strain, "foldchange_", "")
	)

# Plot fold change of PB2 in PA-FS and del-NS1 relative to WT
p_PB2_ratio <- ggplot(
	data_PB2_ratio_long,
	aes(x = factor(strain, levels = c("PA_FS", "del_NS1")), y = FoldChange)) +
	geom_point(
		position = position_jitter(width = 0.1)) +
	labs(
		title = "PB2 in PA-FS and del-NS1 relative to WT",
		x = "", 
		y = "Fold Change") +
	ylim(0,1.1) +
	scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
	theme_classic()
ggsave("TN045/45.8/TN045.8_PB2_WT_ratio_plot.png", plot = p_PB2_ratio, width = 8, height = 6, dpi = 300)
