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

# Statistical testing (test on ddCT; identical to dCT, consistent w/ fold-change framing)
library(rstatix)   # dunnett_test(), levene_test(), t_test()
library(ggpubr)    # stat_pvalue_manual()

ddct_stats <- data_ddct |>
  filter(`Target Name` != "RNA18S") |>
  mutate(condition = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS")))

## DUSP11 (host gene): one-way ANOVA + Tukey HSD (all pairwise, incl. WT vs PA-FS)
dusp11 <- filter(ddct_stats, `Target Name` == "DUSP11")
aov_dusp11 <- aov(ddCT ~ condition, data = dusp11)
print(summary(aov_dusp11))
print(shapiro.test(residuals(aov_dusp11)))     # normality of residuals
print(levene_test(dusp11, ddCT ~ condition))   # equal variances

# Tukey HSD: mock-WT, mock-PA-FS, WT-PA-FS under one family-wise correction
tukey_dusp11 <- tukey_hsd(aov_dusp11)
print(tukey_dusp11)

## WSN_PB2 (viral gene): WT vs PA-FS only; mock is uninfected background. Welch t-test.
pb2 <- ddct_stats |>
  filter(`Target Name` == "WSN_PB2", condition %in% c("WT", "PA-FS")) |>
  mutate(condition = droplevels(condition))
ttest_pb2 <- t_test(pb2, ddCT ~ condition, var.equal = FALSE) |>  # Welch
  add_significance("p")
print(ttest_pb2)

# --- Plot data (per-replicate + means, both targets) ---
plot_data <- data_foldchange |>
  filter(`Target Name` != "RNA18S") |>
  mutate(condition = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS")))

plot_means <- plot_data |>
  group_by(condition, `Target Name`) |>
  summarize(mean_fc = mean(fold_change, na.rm = TRUE),
            sd = sd(fold_change, na.rm = TRUE), .groups = "drop")

# --- DUSP11 plot: linear axis, Dunnett per-bar stars ---
dusp11_means <- filter(plot_means, `Target Name` == "DUSP11")
dusp11_pts   <- filter(plot_data,  `Target Name` == "DUSP11")

# Pairwise brackets (group1 -> group2), staggered so they don't overlap
top_dusp11 <- max(dusp11_means$mean_fc + dusp11_means$sd, na.rm = TRUE)
stars_dusp11 <- tukey_dusp11 |>
  mutate(y.position = top_dusp11 * c(1.1, 1.25, 1.4))

# Colors for plotting — one palette per plot, edit independently
dusp11_colors <- c(
  mock    = "#2A3752",  # navy — echoes the slide header
  WT      = "#C8663E",  # terracotta
  "PA-FS" = "#D9A441"   # amber
)

ggplot(dusp11_means, aes(condition, mean_fc)) +
  geom_col(fill = '#2A3752', alpha = 0.5) +
  geom_jitter(data = dusp11_pts,
              aes(condition, fold_change), color = "#2A3752", width = 0.1) +
  geom_errorbar(aes(ymin = mean_fc - sd, ymax = mean_fc + sd), color = "#2A3752", width = 0.2) +
  stat_pvalue_manual(stars_dusp11, label = "p.adj.signif", tip.length = 0.01) +
  labs(title = "", x = "", y = "Foldchange") +
  scale_x_discrete(labels = c("mock" = "Mock", "WT" = "WT", "PA-FS" = "ΔPA-X")) +
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(), plot.background = element_blank())
ggsave("TN045/45.9/foldchange_DUSP11.png", width = 6, height = 4, dpi = 300, bg = "transparent")

# --- WSN_PB2 plot: log10 axis (viral gene spans ~1e7), WT-vs-PA-FS bracket ---
# SD error bars omitted: symmetric linear SD breaks on a log axis; jitter shows spread.
# mock dropped: uninfected background, not informative for the viral gene.
pb2_means <- plot_means |>
  filter(`Target Name` == "WSN_PB2", condition != "mock") |>
  mutate(condition = droplevels(condition))
pb2_pts <- plot_data |>
  filter(`Target Name` == "WSN_PB2", condition != "mock") |>
  mutate(condition = droplevels(condition))

# y.position is in original data units even under scale_y_log10() -> place above PA-FS
stars_pb2 <- ttest_pb2 |>
  mutate(y.position = max(pb2_pts$fold_change[pb2_pts$condition == "PA-FS"],
                          na.rm = TRUE) * 3)

# PB2 palette (mock dropped) — edit independently of dusp11_colors
pb2_colors <- c(
  WT      = "#2A3752",  # terracotta
  "PA-FS" = "#C8663E"   # amber
)

ggplot(pb2_means, aes(condition, mean_fc)) +
  geom_col(fill = '#2A3752', alpha = 0.5) +
  geom_jitter(data = pb2_pts,
              aes(condition, fold_change, color = "#2A3752"), width = 0.1) +
  stat_pvalue_manual(stars_pb2, x = "group2", label = "p.signif") +  # star over PA-FS
  scale_y_log10() +
  labs(title = "PB2", x = "",
       y = "Fold change (2^-ddCT), log10",
       caption = "Fold change relative to mock; test = WT vs PA-FS, Welch t") +
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank())
ggsave("TN045/45.9/foldchange_WSN_PB2.png", width = 6, height = 4, dpi = 300)
