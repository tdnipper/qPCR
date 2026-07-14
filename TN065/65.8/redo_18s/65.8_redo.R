library(tidyverse)

# Read in data
data <- read_csv("TN065/65.8/redo_18s/TN065.8_redo_filtered.csv")

# Take geo mean of technical reps
data <- data |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  summarize(mean_ct = exp(mean(log(CT), na.rm = TRUE)), .groups = "drop") |>
  pivot_wider(names_from = `Target Name`, values_from = mean_ct)

# calculate dCT
data <- data |>
  mutate(across(c(DUSP11), 
                list(dct = ~ . - `18S`),
                .names = "dct {.col}"))

# pivot longer to get dct columns in one column
data_dct <- data |>
  pivot_longer(cols = starts_with("dct"), names_to = "Target", values_to = "dct_value") |>
  mutate(Target = str_remove(Target, "dct ")) |>
  select(`Sample Name`, Task, Target, dct_value)

# get mean dct of DUSP11 for T0 as control average
control_avg <- data_dct |>
  filter(`Sample Name` == "T0") |>
  group_by(Target) |>
  summarize(control_avg = mean(dct_value, na.rm = TRUE), .groups = "drop")

# join control dct back to data_dct and calculate ddct
data_ddct <- data_dct |>
  left_join(control_avg, by = "Target") |>
  mutate(ddct = dct_value - control_avg)

# get fold change by calculating 2^-ddct
data_foldchange <- data_ddct |>
  mutate(fold_change = 2^(-ddct))

# save results to csv
write_csv(
  data_foldchange |>
    mutate(across(c(dct_value, control_avg, ddct, fold_change),
    ~ signif(.x, 2)))
  ,"TN065/65.8/redo_18s/foldchange_results.csv"
)

# --- Statistical testing: one-way ANOVA + Dunnett vs T0 ---
library(rstatix)   # dunnett_test(), levene_test()

# Test on dCT (log2 scale), not fold change, which is log-normal / heteroscedastic.
# T0 first in levels = Dunnett reference/control.
dct_stats <- data_dct |>
  filter(Target == "DUSP11") |>
  mutate(timepoint = factor(`Sample Name`, levels = c("T0", "T8", "T24", "T48")))

aov_dusp11 <- aov(dct_value ~ timepoint, data = dct_stats)
print(summary(aov_dusp11))

# Assumption checks (n=6/group -> ANOVA fairly robust)
print(shapiro.test(residuals(aov_dusp11)))            # normality of residuals
print(levene_test(dct_stats, dct_value ~ timepoint))  # equal variances

# Dunnett: each timepoint vs T0 (single-step FWER correction)
dunnett <- dunnett_test(dct_stats, dct_value ~ timepoint, ref.group = "T0")
print(dunnett)

# data to plot foldchange
plot_data <- data_foldchange |>
  select(`Sample Name`, Task, Target, fold_change)

# order x axis by timepoint
plot_data <- plot_data |>
  mutate(`Sample Name` = factor(`Sample Name`, levels = c("T0", "T8", "T24", "T48")))

# take means of fold change for each timepoint and target
plot_data_means <- plot_data |>
  group_by(`Sample Name`, Target) |>
  summarize(
    mean_fold_change = mean(fold_change, na.rm = TRUE),
    sd = sd(fold_change, na.rm = TRUE),
    .groups = "drop")

# place the Dunnett stars above each error bar (stars come from dunnett_test's p.adj.signif)
library(ggpubr)    # stat_pvalue_manual()

stat_test <- dunnett |>
  mutate(`Sample Name` = factor(group2, levels = c("T0", "T8", "T24", "T48"))) |>
  left_join(plot_data_means, by = "Sample Name") |>
  mutate(y.position = mean_fold_change + sd + 0.05)

# plot fold change
ggplot(plot_data_means, aes(x = `Sample Name`, y = mean_fold_change)) +
  geom_col(aes(fill = `Sample Name`), alpha = 0.5) +
  geom_jitter(data = plot_data, aes(x = `Sample Name`, y = fold_change, color = `Sample Name`), width = 0.2) +
  geom_errorbar(aes(ymin = mean_fold_change - sd,
                ymax = mean_fold_change + sd), width = 0.2) +
  stat_pvalue_manual(stat_test, x = "group2", label = "p.adj.signif") +
  labs(title = "DUSP11 during infection", x = "Hours post-infection", y = "Fold change (2^-ddCT)") +
  scale_y_continuous(limits = c(0, 1.4), breaks = seq(0, 1.25, by = 0.25)) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave("TN065/65.8/redo_18s/foldchange_plot.png", width = 6, height = 4, dpi = 300)