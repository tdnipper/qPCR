library(tidyverse)

data <- read_csv("TN071/71.1/TN071.1_filtered.csv")

data <- data |>
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
		.cols = targets, 
		.fn = ~ . - TBP1, 
		.names = "dCT_{col}")) |>
  ungroup()

# Calculate mock mean dCT of DUSP11 and PB2 and store as named vector
mock_values <- data_dct |>
  filter(`Sample Name` == "mock") |>
  summarize(
	across(
		.cols = starts_with("dCT_"), 
		mean)) |>
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

# Plot fold change for DUSP11 and PB2

p_DUSP11 <- ggplot(data_foldchange, aes(x = `Sample Name`, y = foldchange_DUSP11)) +
  geom_point(position = position_jitter(width = 0.1)) +
  labs(
	title = "DUSP11 expression during infection",
	x = "",
	y = "Fold Change") +
  theme_classic()
ggsave("TN071/71.1/TN071.1_DUSP11_foldchange_plot.png", plot = p_DUSP11, width = 8, height = 6, dpi = 300)

p_PB2 <- ggplot(data_foldchange, aes(x = `Sample Name`, y = foldchange_WSN_PB2)) +
  geom_point(position = position_jitter(width = 0.1)) +
  labs(
	title = "PB2 expression during infection",
	x = "",
	y = "Fold Change") +
  scale_y_log10() +
  theme_classic()
ggsave("TN071/71.1/TN071.1_PB2_foldchange_plot.png", plot = p_PB2, width = 8, height = 6, dpi = 300)

# Plot amplification data for housekeeping gene for each group
data_amp <- read_csv("TN071/71.1/TN071.1_amp.csv") |>
  filter(`Target Name` == "TBP1") |>
  group_by(
	`Sample Name`,
	Cycle,
	Task) |>
  summarize(
	meanRn = mean(`Delta Rn`), 
	sdRN = sd(`Delta Rn`),
	.groups = "drop")

p_amp <- ggplot(
	data_amp, 
	aes(
		x = Cycle, 
		y = meanRn, 
		color = `Sample Name`, 
		group = interaction(
			`Sample Name`, 
			Task)
	)
  ) +
  geom_line() +
  geom_ribbon(
	aes(
	ymin = meanRn - sdRN,
	ymax = meanRn + sdRN,
	fill = `Sample Name`
	),
	alpha = 0.2,
  ) +
  labs(
	title = "TBP Amplification",
	x = "",
	y = "CT") +
  theme_classic()

ggsave("TN071/71.1/TN071.1_TBP1_amplification_plot.png", plot = p_amp, width = 8, height = 6, dpi = 300)