library(tidyverse)
library(rstatix)

df <- read_csv("TN072/72.4/TN072.4_filtered.csv", name_repair = "universal") 

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

# Calculate dCT for both targets (normalized to RNA18S1)
df_dct <- df |>
  pivot_wider(names_from = Target.Name, values_from = mean_ct) |>
  mutate(dct_DUSP11 = DUSP11 - RNA18S1,
         dct_ISG15 = ISG15 - RNA18S1)

# ddCT for both targets, using mean mock dCT as control
control_avg <- df_dct |>
  filter(Sample.Name == "mock") |>
  summarize(across(c(dct_DUSP11, dct_ISG15), mean))

df_ddct <- df_dct |>
  mutate(ddct_DUSP11 = dct_DUSP11 - control_avg$dct_DUSP11,
         ddct_ISG15 = dct_ISG15 - control_avg$dct_ISG15,
         foldchange_DUSP11 = 2^(-ddct_DUSP11),
         foldchange_ISG15 = 2^(-ddct_ISG15))

print(df_ddct, n = Inf)


# Long format: one row per sample x replicate x gene
df_long <- df_ddct |>
  select(Sample.Name, Task, ddct_DUSP11, ddct_ISG15) |>
  pivot_longer(starts_with("ddct_"), names_to = "gene", names_prefix = "ddct_", values_to = "ddct") |>
  mutate(Sample.Name = factor(Sample.Name, levels = c("mock", "IFN", "polyI:C", "infect")))

# Stats on ddCT (log2 scale), per gene: ANOVA, assumption checks, Dunnett vs mock
cat("\n== ANOVA ==\n")
df_long |> group_by(gene) |> anova_test(ddct ~ Sample.Name) |> print()

cat("\n== Shapiro on ANOVA residuals ==\n")
df_long |>
  group_by(gene) |>
  group_modify(~ broom::tidy(shapiro.test(residuals(aov(ddct ~ Sample.Name, data = .x))))) |>
  print()

cat("\n== Levene ==\n")
df_long |> group_by(gene) |> levene_test(ddct ~ Sample.Name) |> print()

set.seed(42)
stat_test <- df_long |>
  group_by(gene) |>
  dunnett_test(ddct ~ Sample.Name, ref.group = "mock") |>
  add_significance("p.adj")
cat("\n== Dunnett vs mock ==\n")
print(stat_test, width = Inf)

# Plot data: mean and SD on ddCT, back-transformed to fold change (asymmetric bounds)
df_plot <- df_long |>
  summarize(mean_ddct = mean(ddct), sd_ddct = sd(ddct), n = n(), .by = c(gene, Sample.Name)) |>
  mutate(mean_fc = 2^(-mean_ddct),
         fc_lower = 2^(-(mean_ddct + sd_ddct)),
         fc_upper = 2^(-(mean_ddct - sd_ddct)))
print(df_plot, n = Inf, width = Inf)


# Bar chart: fold change vs mock, stars from Dunnett on ddCT
stat_plot <- stat_test |>
  left_join(df_plot |> select(gene, group2 = Sample.Name, fc_upper), by = c("gene", "group2")) |>
  mutate(y.position = fc_upper * 1.05)

p <- ggplot(df_plot, aes(x = Sample.Name, y = mean_fc)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = fc_lower, ymax = fc_upper), width = 0.2) +
  ggpubr::stat_pvalue_manual(stat_plot, x = "group2", label = "p.adj.signif", hide.ns = FALSE) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "DUSP11 and ISG15 expression", x = "Sample Name", y = "Fold change vs mock") +
  theme_minimal()
ggsave("TN072/72.4/TN072.4_foldchange.png", plot = p, width = 8, height = 4, dpi = 300)
