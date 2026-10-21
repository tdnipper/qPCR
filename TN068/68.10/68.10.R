library(tidyverse)
library(rstatix)

df <- read_csv("TN068/68.10/TN068.10_filtered.csv", name_repair = "universal") 

# Trim outlier furthest from mean in SD > threshold, but only if the remaining wells agree (SD <= threshold)
threshold <- 0.5

is_outlier <- function(ct) {
  worst <- which.max(abs(ct - median(ct)))
  trim <- length(ct) >= 3 && sd(ct) > threshold && sd(ct[-worst]) <= threshold
  seq_along(ct) == worst & trim
}

df_flagged <- df |>
  filter(!is.na(CT)) |>
  mutate(outlier = is_outlier(CT), .by = c(Sample.Name, Target.Name, Task))

cat("\n== technical replicates dropped ==\n")
df_flagged |> filter(outlier) |> print(n = Inf)

# Collapse technical reps on the survivors
df <- df_flagged |>
  filter(!outlier) |>
  summarize(mean_ct = mean(CT), sd_ct = sd(CT), .by = c(Sample.Name, Target.Name, Task))

cat("\n== unresolved groups (kept untrimmed — review/re-run) ==\n")
df |> filter(sd_ct > threshold) |> print(n = Inf)

df <- df |> select(-sd_ct)

# Test: mock samples only, Dox vs noDox (dCT, per cell line and target)
df_mock <- df |>
  filter(str_detect(Sample.Name, "mock")) |>
  pivot_wider(names_from = Target.Name, values_from = mean_ct) |>
  mutate(
    Cell.Line = str_extract(Sample.Name, "^[^_]+"),
    Dox = if_else(str_detect(Sample.Name, "noDox"), "noDox", "Dox") |> factor(levels = c("noDox", "Dox")),
    across(c(ISG15, MX1), ~ .x - RNA18S1, .names = "dct_{col}")
  )

# ddCT of Dox vs the mean noDox mock dCT of the same cell line
df_mock <- df_mock |>
  mutate(
    across(starts_with("dct_"), ~ .x - mean(.x[Dox == "noDox"], na.rm = TRUE), .names = "d{.col}"),
    .by = Cell.Line
  )
print(df_mock |> select(Sample.Name, Task, starts_with("d")), n = Inf)

mock_test <- df_mock |>
  pivot_longer(starts_with("dct_"), names_to = "Target", values_to = "dct") |>
  group_by(Cell.Line, Target) |>
  t_test(dct ~ Dox) |>
  add_significance()
cat("\n== mock only: Dox vs noDox (t-test on dCT) ==\n")
print(mock_test)

# Fold change (2^-ddCT, Dox vs noDox mock) - mean and SD across replicates
fc_summary <- df_mock |>
  pivot_longer(starts_with("ddct_"), names_to = "Target", values_to = "ddct", names_prefix = "ddct_") |>
  mutate(fc = 2^-ddct) |>
  summarize(mean_fc = mean(fc), sd_fc = sd(fc), .by = c(Cell.Line, Target, Dox))

p_fc <- fc_summary |>
  ggplot(aes(x = Target, y = mean_fc, fill = Dox)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc),
                position = position_dodge(width = 0.8), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~Cell.Line) +
  labs(x = NULL, y = "Fold change vs noDox mock (2^-ddCT)") +
  theme_minimal()
ggsave("TN068/68.10/68.10_mock_dox_foldchange_plot.png", p_fc, dpi = 300, height = 4, width = 6)

# Comparison 2: infect vs mock, within each cell line and Dox condition
df_inf <- df |>
  pivot_wider(names_from = Target.Name, values_from = mean_ct) |>
  mutate(
    Cell.Line = str_extract(Sample.Name, "^[^_]+"),
    Dox = if_else(str_detect(Sample.Name, "noDox"), "noDox", "Dox") |> factor(levels = c("noDox", "Dox")),
    Infection = if_else(str_detect(Sample.Name, "mock"), "mock", "infect") |> factor(levels = c("mock", "infect")),
    across(c(ISG15, MX1), ~ .x - RNA18S1, .names = "dct_{col}")
  ) |>
  mutate(
    across(starts_with("dct_"), ~ .x - mean(.x[Infection == "mock"], na.rm = TRUE), .names = "d{.col}"),
    .by = c(Cell.Line, Dox)
  )

inf_test <- df_inf |>
  pivot_longer(starts_with("dct_"), names_to = "Target", values_to = "dct", names_prefix = "dct_") |>
  group_by(Cell.Line, Dox, Target) |>
  t_test(dct ~ Infection) |>
  add_significance()
cat("\n== infect vs mock (t-test on dCT) ==\n")
print(inf_test)

inf_fc <- df_inf |>
  filter(Infection == "infect") |>
  pivot_longer(starts_with("ddct_"), names_to = "Target", values_to = "ddct", names_prefix = "ddct_") |>
  mutate(fc = 2^-ddct) |>
  summarize(mean_fc = mean(fc), sd_fc = sd(fc), .by = c(Cell.Line, Dox, Target))

p_inf <- inf_fc |>
  ggplot(aes(x = Target, y = mean_fc, fill = Dox)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc),
                position = position_dodge(width = 0.8), width = 0.2) +
  facet_wrap(~Cell.Line) +
  labs(x = NULL, y = "Fold change vs mock (2^-ddCT)") +
  theme_minimal()
ggsave("TN068/68.10/68.10_infect_vs_mock_foldchange_plot.png", p_inf, dpi = 300, height = 4, width = 6)
