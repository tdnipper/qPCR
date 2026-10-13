library(tidyverse)
library(rstatix)

df <- read_csv("TN065/65.9/TN065.9_filtered.csv", name_repair = "universal")

# Technical-replicate outliers: if a group's SD > threshold, drop a well only when
# it is the *only* one whose removal brings the rest within threshold
threshold <- 0.5

is_outlier <- function(ct) {
  fixes <- map_lgl(seq_along(ct), \(i) isTRUE(sd(ct[-i]) <= threshold))
  sd(ct) > threshold & sum(fixes) == 1 & fixes
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

df <- df |>
  select(-sd_ct) |>
  pivot_wider(names_from = Target.Name, values_from = mean_ct)

# Calculate dCT = mean_ct - mean_ct(mock) for each target, using mock as control
df <- df |>
  mutate(across(c(DUSP11, WSN_PB2), 
                 ~ . - `RNA18S1`,
                .names = "dct_{.col}"))

df_dct <- df |>
  pivot_longer(cols = starts_with("dct_"), names_to = "Target", values_to = "dct_value") |>
  mutate(Target = str_remove(Target, "dct_")) |>
  select(Sample.Name, Task, Target, dct_value)

# Get mean of dCT at mock as control average
control_avg <- df_dct |>
  filter(Sample.Name == "mock") |>
  group_by(Target) |>
  summarize(control_avg = mean(dct_value, na.rm = TRUE), .groups = "drop")

df_ddct <- df_dct |>
  left_join(control_avg, by = "Target") |>
  mutate(ddct = dct_value - control_avg) |>
  mutate(fold_change = 2^(-ddct))

write_csv(df_ddct, "TN065/65.9/TN065.9_ddct.csv")

# Stats: test on ddCT
timepoints <- c("mock", "T0", "T1", "T2", "T4", "T6", "T8", "T24")

ddct_stats <- df_ddct |>
  filter(Target == "DUSP11", !is.na(ddct)) |>
  mutate(Sample.Name = factor(Sample.Name, levels = timepoints))

aov_dusp11 <- aov(ddct ~ Sample.Name, data = ddct_stats)
print(summary(aov_dusp11))

# Assumption checks
print(shapiro.test(residuals(aov_dusp11)))          # normality of residuals
print(levene_test(ddct_stats, ddct ~ Sample.Name))  # equal variances

# Dunnett: each timepoint vs mock
set.seed(42)
dunnett <- dunnett_test(ddct_stats, ddct ~ Sample.Name, ref.group = "mock")
print(dunnett)

df_plot <- df_ddct |>
  filter(Target == "DUSP11") |>
  group_by(Sample.Name, Target) |>
  summarize(
    mean_fold_change = mean(fold_change, na.rm = TRUE),
    sd_fold_change = sd(fold_change, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(Sample.Name = factor(Sample.Name, levels = timepoints))

# place the Dunnett stars above each error bar (stars come from dunnett_test's p.adj.signif)
library(ggpubr)    # stat_pvalue_manual()

stat_test <- dunnett |>
  mutate(Sample.Name = factor(group2, levels = timepoints)) |>
  left_join(df_plot, by = "Sample.Name") |>
  mutate(y.position = mean_fold_change + coalesce(sd_fold_change, 0) + 0.05 * max(df_plot$mean_fold_change))

p <- ggplot(
  df_plot,
  aes(x = Sample.Name, y = mean_fold_change)
) +
  geom_col(fill = "#2A3752", alpha = 0.8) +
  geom_errorbar(
    color = "#2A3752",
    aes(
      ymin = mean_fold_change - sd_fold_change,
      ymax = mean_fold_change + sd_fold_change
    ),
    width = 0.2,
    na.rm = TRUE
  ) +
  stat_pvalue_manual(stat_test, x = "group2", label = "p.adj.signif", hide.ns = TRUE) +
  labs(
    title = "DUSP11 during sychronized single cycle infection",
    y = "Foldchange",
    x = "",
    caption = "MOI 5"
    ) +
  theme_minimal() +
  theme(panel.grid = element_blank())
print(p)
ggsave(
  filename = "TN065/65.9/65.9_foldchange.png",
  plot = p,
  dpi = 300,
  height = 4,
  width = 6
)
