library(tidyverse)

# Read in data
data <- read_csv("TN045/45.9/TN045.9_filtered.csv")

# Take geometric mean of technical reps
data <- data |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  summarize(ct = exp(mean(log(CT), na.rm = TRUE)), .groups = "drop")

# Get non housekeeping targets
targets <- data |>
  filter(`Target Name` != "RNA18S") |>
  pull(`Target Name`) |>
  unique()

# Calculate dCT
data_dct <- data |>
  group_by(`Sample Name`, Task) |>
  pivot_wider(names_from = `Target Name`, values_from = ct) |>
  mutate(
    across(
      .cols = all_of(targets),
      .fn = ~ . - RNA18S,
      .names = "dCT_{col}"
    )
  ) |>
  ungroup()

data_dct <- data_dct |>
  pivot_longer(cols = starts_with("dCT_"), names_to = "Target Name", values_to = "dCT") |>
  mutate(`Target Name` = str_remove(`Target Name`, "dCT_")) |>
  select(`Sample Name`, Task, `Target Name`, dCT)

# Calculate control average from mock
control_avg <- data_dct |>
  filter(`Sample Name` == "mock") |>
  group_by(`Target Name`) |>
  summarize(control_avg = mean(dCT, na.rm = TRUE), .groups = "drop")

# Calculate ddCT
data_ddct <- data_dct |>
  left_join(control_avg, by = "Target Name") |>
  mutate(ddCT = dCT - control_avg)

# calculate fold change as 2^(-ddCT)
data_foldchange <- data_ddct |>
  mutate(fold_change = 2^(-ddCT))

# save results to csv
write_csv(
  data_foldchange |>
    mutate(across(c(dCT, control_avg, ddCT, fold_change),
    ~ signif(.x, 2)))
  ,"TN045/45.9/foldchange_results.csv"
)

# Statistical testing: one-way ANOVA + Dunnett vs mock
library(rstatix)   # dunnett_test(), levene_test()
set.seed(123)  # for reproducibility

# test on ddct
ddct_stats <- data_ddct |>
  filter(`Target Name` != "RNA18S") |>
  mutate(timepoint = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS")))

