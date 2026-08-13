library(tidyverse)

# Read in data
data <- read_csv("TN045/45.10/TN045.10_filtered.csv")

# Flag outlier technical reps: for groups with sd > threshold, mark the single
# replicate furthest from the median (robust to the outlier itself). Same rule as
# TN074/74.1, with two guards: never trim a group that is already a duplicate, and
# never drop more than one well per group.
threshold <- 0.5

data_flagged <- data |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  mutate(
    med        = median(CT, na.rm = TRUE),
    dist       = abs(CT - med),
    grp_sd     = sd(CT, na.rm = TRUE),
    is_outlier = n() >= 3 & grp_sd > threshold & dist == max(dist) & CT != med
  ) |>
  ungroup()

# inspect which wells are dropped before trusting downstream numbers
cat("\n== technical replicates dropped (sd >", threshold, "CT) ==\n")
data_flagged |>
  filter(is_outlier) |>
  select(`Sample Name`, `Target Name`, Task, CT, med, grp_sd) |>
  print(n = Inf)

# Collapse technical reps on the survivors. CT is already log2(quantity), so the
# arithmetic mean of CT *is* the geometric mean in copy-number space.
data <- data_flagged |>
  filter(!is_outlier) |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  summarize(
    n_reps = n(),
    ct     = mean(CT, na.rm = TRUE),
    sd_ct  = sd(CT, na.rm = TRUE),
    .groups = "drop"
  )

# groups still over threshold after one removal are unreliable — review/re-run
cat("\n== groups still over threshold after removal (unreliable) ==\n")
data |> filter(sd_ct > threshold) |> print(n = Inf)

data <- data |> select(-n_reps, -sd_ct)

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

# Replicate-level outliers. Distinct from the technical-rep step above: a whole
# replicate can be internally consistent (tight triplicate) yet sit far from its
# condition. No technical-rep rule can catch these — mock Task 5 reads 21.17 /
# 22.27 / 21.19, so every choice of which well to drop still leaves it ~1 CT low.
# Flag on robust MAD z of dCT within condition; mad() already applies the 1.4826
# scaling, so do NOT multiply by it again.
z_threshold <- 3.5

data_dct <- data_dct |>
  group_by(`Sample Name`, `Target Name`) |>
  mutate(
    mad_dCT = mad(dCT, na.rm = TRUE),
    mad_z   = if_else(mad_dCT > 0,
                      abs(dCT - median(dCT, na.rm = TRUE)) / mad_dCT,
                      NA_real_),
    is_rep_outlier = !is.na(mad_z) & mad_z > z_threshold
  ) |>
  ungroup()

cat("\n== replicate-level outliers excluded (MAD z >", z_threshold, ") ==\n")
data_dct |>
  filter(is_rep_outlier) |>
  select(`Sample Name`, `Target Name`, Task, dCT, mad_z) |>
  print(n = Inf)

data_dct <- data_dct |>
  filter(!is_rep_outlier) |>
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
  mutate(sample = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS")))

# DUSP11 one way ANOVA + Tukey HSD (all pairwise, incl. WT vs PA-FS)
dusp11 <- filter(ddct_stats, `Target Name` == "DUSP11")
aov_dusp11 <- aov(ddCT ~ sample, data = dusp11)
print(summary(aov_dusp11))
print(shapiro.test(aov_dusp11$residuals))
print(levene_test(dusp11, ddCT ~ sample))
print(tukey_hsd(dusp11, ddCT ~ sample))

welch_dusp11 <- welch_anova_test(dusp11, ddCT ~ sample)
gh_dusp11    <- games_howell_test(dusp11, ddCT ~ sample)
print(welch_dusp11)
print(gh_dusp11)

# plot data
plot_data <- data_foldchange |>
  filter(`Target Name` == "DUSP11") |>
  mutate(`Sample Name` = factor(`Sample Name`, levels = c("mock", "WT", "PA-FS")))

plot_means <- plot_data |>
  group_by(`Sample Name`) |>
  summarize(mean_ddCT = mean(ddCT, na.rm = TRUE),
            sd_ddCT = sd(ddCT, na.rm = TRUE), .groups = "drop") |>
  mutate(mean_fc = 2^(-mean_ddCT),
         ymin = 2^(-(mean_ddCT + sd_ddCT)),
         ymax = 2^(-(mean_ddCT - sd_ddCT)))

sample_levels <- c("mock", "WT", "PA-FS")
top_dusp11 <- max(c(plot_data$fold_change, plot_means$ymax), na.rm = TRUE)

stars_dusp11 <- gh_dusp11 |>
  arrange(match(group1, sample_levels), match(group2, sample_levels)) |>
  mutate(y.position = top_dusp11 * (1.10 + 0.15 * (row_number() - 1)))

# DUSP11 plot
ggplot() +
  geom_col(data = plot_means, aes(x = `Sample Name`, y = mean_fc), fill = "#2A3752", alpha = 0.5) +
  geom_jitter(data = plot_data, aes(x = `Sample Name`, y = fold_change), color = "#2A3752", width = 0.1) +
  geom_errorbar(data = plot_means, aes(x = `Sample Name`, ymin = ymin, ymax = ymax), color = "#2A3752", width = 0.2) +
  stat_pvalue_manual(stars_dusp11, label = "p.adj.signif", tip.length = 0.01) +
  labs(title = "",
       x = "", y = "Foldchange") +
  theme_minimal()

ggsave("TN045/45.10/DUSP11_fold_change_plot.png", width = 6, height = 4, dpi = 300)