library (tidyverse)

# read in data
data <- read_csv("TN074/74.1/TN074.1.2_filtered.csv")

# flag outlier technical reps: for groups with sd > threshold, mark the single
# replicate furthest from the median (robust to the outlier itself)
threshold <- 0.5

data_flagged <- data |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  mutate(
    med        = median(CT, na.rm = TRUE),
    dist       = abs(CT - med),
    grp_sd     = sd(CT, na.rm = TRUE),
    is_outlier = dist == max(dist) & CT != med & grp_sd > threshold
  ) |>
  ungroup()

# inspect which wells would be dropped before committing
data_flagged |> filter(is_outlier) |> print(n = Inf)

# drop the flagged wells, then recompute means and sd on the survivors
data_mean <- data_flagged |>
  filter(!is_outlier) |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  summarize(
    n_reps  = n(),
    mean_Ct = exp(mean(log(CT), na.rm = TRUE)),
    sd      = sd(CT, na.rm = TRUE),
    .groups = "drop"
  )

# groups still above threshold after cleanup are unreliable — review/re-run
data_mean |> filter(sd > threshold) |> print(n = Inf)

# Pivot wider and calculate dCT
data_dCT <- data_mean |>
  select(-sd, -n_reps) |>
  pivot_wider(
    names_from  = `Target Name`,
    values_from = mean_Ct
  ) |>
  mutate(
    dCT = `DUSP11` - `RNA18S1`
  )

# calculate mean control dct
control_avg <- data_dCT |>
  filter(`Sample Name` == "mock") |>
  summarize(control_avg = mean(dCT, na.rm = TRUE), .groups = "drop")

# join control dct back to data_dCT and calculate ddCT
data_ddCT <- data_dCT |>
  cross_join(control_avg) |>
  mutate(ddCT = dCT - control_avg) |>
  select(`Sample Name`, Task, ddCT)

# get fold change by calculating 2^-ddCT
data_fold_change <- data_ddCT |>
  mutate(fold_change = 2^(-ddCT))

# write out the fold change data
write_csv(
  mutate(
    data_fold_change,
    across(
      .cols = c(ddCT, fold_change),
      .fns  = ~ signif(.x, digits = 2)
    )
  ),
  "TN074/74.1/TN074.1.2_fold_change.csv"
)

# get means to plot
data_plot <- data_fold_change |>
  group_by(`Sample Name`) |>
  summarize(
    mean_fold_change = mean(fold_change, na.rm = TRUE),
    sd_fold_change   = sd(fold_change, na.rm = TRUE),
    .groups = "drop"
  )

# plot fold change data
ggplot(data_plot, aes(x = factor(`Sample Name`, levels = c("mock", "infected")), y = mean_fold_change)) +
  geom_col(alpha = 0.7) +
  geom_errorbar(
    aes(ymin = mean_fold_change - sd_fold_change, ymax = mean_fold_change + sd_fold_change),
    width = 0.2
  ) +
  geom_jitter(data = data_fold_change, aes(x = `Sample Name`, y = fold_change, color = factor(`Task`)), width = 0.01) +
  labs(
    title    = "DUSP11 Infection 2",
    x        = "",
    y        = "Fold Change",
    color = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 1.25), breaks = seq(0, 1.25, by = 0.25))
ggsave("TN074/74.1/TN074.1.2_foldchange_plot.png", width = 6, height = 4, units = "in", dpi = 300)