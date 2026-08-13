library(tidyverse)

# Read in data
data <- read_csv("TN045/45.10/TN045.10_filtered.csv")

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
  ,"TN045/45.10/foldchange_results.csv"
)

# Statistical testing (test on ddCT; identical to dCT, consistent w/ fold-change framing)
library(rstatix)   # dunnett_test(), levene_test(), t_test()
library(ggpubr)    # stat_pvalue_manual()

ddct_stats <- data_ddct |>
  filter(`Target Name` != "RNA18S") |>
  mutate(condition = factor(`Target Name`, levels = c("mock", "WT", "PA-FS")))

# DUSP11 one way ANOVA + Tukey HSD (all pairwise, incl. WT vs PA-FS)
dusp11 <- filter(ddct_stats, `Target Name` == "DUSP11")
aov_dusp11 <- aov(ddCT ~ `Sample Name`, data = dusp11)
print(summary(aov_dusp11))
print(shapiro.test(aov_dusp11$residuals))
print(levene_test(dusp11, ddCT ~ `Sample Name`))

# plot data
plot_data <- data_foldchange |>
  filter(`Target Name` == "DUSP11") |>
  mutate(`Sample Name` = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS")))

plot_means <- plot_data |>
  group_by(`Sample Name`) |>
  summarize(mean_fc = mean(fold_change, na.rm = TRUE),
            sd_fc = sd(fold_change, na.rm = TRUE), .groups = "drop")

# Colors for plotting — one palette per plot, edit independently
dusp11_colors <- c(
  mock    = "#2A3752",  # navy — echoes the slide header
  WT      = "#C8663E",  # terracotta
  "PA-FS" = "#D9A441"   # amber
)

# DUSP11 plot
ggplot() +
  geom_col(data = plot_means, aes(x = `Sample Name`, y = mean_fc), fill = "#2A3752", alpha = 0.5) +
  geom_jitter(data = plot_data, aes(x = `Sample Name`, y = fold_change), color = "#2A3752", width = 0.1) +
  geom_errorbar(data = plot_means, aes(x = `Sample Name`, ymin = mean_fc - sd_fc, ymax = mean_fc + sd_fc), color = "#2A3752", width = 0.2) +
  labs(title = "DUSP11 Fold Change", x = "Sample Name", y = "Foldchange (2^-ddCT)") +
  scale_fill_manual(values = dusp11_colors) +
  theme_minimal()

ggsave("TN045/45.10/DUSP11_fold_change_plot.png", width = 6, height = 4, dpi = 300)