library(tidyverse)

# read in data
data <- read_csv("TN074/74.1/TN074.1.1_filtered.csv")

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
data_foldchange <- data_ddCT |>
  mutate(fold_change = 2^(-ddCT)) |>
  mutate(`Sample Name` = factor(`Sample Name`, levels = c("mock", "infected")))

# save results to csv
write_csv(
  data_foldchange |>
    mutate(across(.cols = c(ddCT, fold_change),
    .fns = ~ signif(.x, 2)))
  ,"TN074/74.1/1.1_foldchange_results.csv"
)

# get means to plot
data_plot <- data_foldchange |>
  group_by(`Sample Name`) |>
  summarize(
    mean_fc = mean(fold_change, na.rm = TRUE),
    sd_fc   = sd(fold_change, na.rm = TRUE),
    .groups = "drop"
  )

# plot means and initial fold change values
ggplot(data_plot, aes(x = `Sample Name`, y = mean_fc)) +
  geom_col() +
  geom_jitter(data = data_foldchange, aes(x = `Sample Name`, y = fold_change, color = factor(Task)), width = 0.1) +
  geom_errorbar(aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc), width = 0.2) +
  labs(x = "Sample Name", y = "Fold Change", title = "DUSP11 during infection", caption = "Infection 1, reps 1, 2", col = "Replicate") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  theme(legend.position = "top") +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.25))
ggsave("TN074/74.1/1.1_foldchange_plot.png", width = 8, height = 6, dpi = 300)