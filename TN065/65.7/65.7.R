library(tidyverse)

data <- read_csv("TN065/65.7/TN065.7_filtered.csv")

data <- data |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  summarize(mean_ct = exp(mean(log(CT), na.rm = TRUE)), .groups = "drop") |>
  pivot_wider(names_from = `Target Name`, values_from = mean_ct)

# get geometric mean of PUM1 and TBP for combined housekeeping control
data <- data |>
  rowwise() |>
  mutate(geo_mean = exp(mean(log(c_across(c(PUM1, TBP)))))) |>
  ungroup()

# create dct columns: subtract specified cols geo mean of PUM1 and TBP
cols_to_subtract <- c("DUSP11", "WSN_PB2")
for (col in cols_to_subtract) {
  newname <- paste0("dct ", col)
  data[[newname]] <- data[[col]] - data$geo_mean
}

# pivot longer to get dct columns in one column
data_dct <- data |>
  pivot_longer(cols = starts_with("dct"), names_to = "Target", values_to = "dct_value") |>
  mutate(Target = str_remove(Target, "dct ")) |>
  select(`Sample Name`, Task, Target, dct_value)

# get mean dct of DUSP11 and WSN_PB2 for T0 as control
control_dct <- data_dct |>
  filter(`Sample Name` == "T0") |>
  group_by(Target) |>
  summarize(control_dct = mean(dct_value, na.rm = TRUE), .groups = "drop")

# join control dct back to data_dct and calculate ddct
data_ddct <- data_dct |>
  left_join(control_dct, by = "Target") |>
  mutate(ddct = dct_value - control_dct)

# get fold change by calculating 2^-ddct
data_foldchange <- data_ddct |>
  mutate(fold_change = 2^(-ddct))

# save results to csv
write_csv(data_foldchange, "TN065/65.7/dct_results.csv")

# data to plot foldchange
plot_data <- data_foldchange |>
  select(`Sample Name`, Task, Target, fold_change)

# order x axis by sample name
plot_data$`Sample Name` <- factor(plot_data$`Sample Name`, levels = c("T0", "T8", "T24", "T48"))
plot_data$Task <- as.factor(plot_data$Task)
# take means of foldchange for geom_col later
means <- plot_data |>
  group_by(`Sample Name`, Target) |>
  summarize(mean = mean(fold_change, na.rm = TRUE),
            sd = sd(fold_change, na.rm = TRUE),
            .groups = "drop")

# plot fold change
p_dusp11 <- plot_data |>
  filter(Target == "DUSP11") |>
  ggplot(aes(x = `Sample Name`, y = fold_change, color = `Sample Name`)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_col(data = means |> filter(Target == "DUSP11"), aes(x = `Sample Name`, y = mean, fill = `Sample Name`), alpha = 0.5, inherit.aes = FALSE) +
    geom_errorbar(data = means |> filter(Target == "DUSP11"), aes(x = `Sample Name`, ymin = mean - sd, ymax = mean + sd, color = `Sample Name`), width = 0.2, inherit.aes = FALSE) +
    labs(title = "DUSP11 During Infection",
        x = "",
        y = "Fold Change (2^-ddCT)",
        color = NULL,
        fill = NULL
        ) +
    theme_classic()
ggsave("TN065/65.7/dusp11_plot.png", plot = p_dusp11, width = 10, height = 6, dpi = 300)

# remove outliers (Task 4 and Task 6) from DUSP11 plot
p_dusp11_no_outliers <- plot_data |>
  filter(Task != 4, Task != 6, Target == "DUSP11") |>
  ggplot(aes(x = `Sample Name`, y = fold_change, color = Task)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "DUSP11 During Infection",
        subtitle = "Outliers (4 and Task 6) Removed",
        x = "",
        y = "Fold Change (2^-ddCT)"
        ) +
    theme_classic()
ggsave("TN065/65.7/dusp11_no_outliers_plot.png", plot = p_dusp11_no_outliers, width = 10, height = 6, dpi = 300)

p_pb2 <- plot_data |>
  filter(Target == "WSN_PB2") |>
  ggplot(aes(x = `Sample Name`, y = fold_change, color = Task)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "PB2 During Infection",
        x = "",
        y = "Fold Change (2^-ddCT)"
        ) +
    scale_y_log10() +
    theme_classic()
ggsave("TN065/65.7/pb2_plot.png", plot = p_pb2, width = 10, height = 6, dpi = 300)