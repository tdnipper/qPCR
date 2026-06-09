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
data_plot <- data_foldchange |>
  pivot_longer(cols = starts_with("foldchange_"), names_to = "Target", values_to = "FoldChange") |>
  mutate(Target = str_replace(Target, "foldchange_", ""))

p <- ggplot(
	data_plot, 
	aes(
		x = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS", "del-NS1")),
		y = FoldChange)
		) +
	geom_point(
	stat = "identity", 
	position = position_jitter(width = 0.1)
	) +
  labs(
		x = "", 
		y = "Fold Change") +
  theme_classic() + 
  facet_wrap(vars(Target), scales = "free_y")

ggsave("TN045/45.8/TN045.8_foldchange_plot.png", plot = p, width = 8, height = 6, dpi = 300)
