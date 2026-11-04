library(tidyverse)
library(rstatix)

df <- read_csv("TN072/72.5/TN072.5_filtered.csv", name_repair = "universal")

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
  summarize(mean_ct = mean(CT), sd_ct = sd(CT), n_wells = n(), .by = c(Sample.Name, Target.Name, Task))

cat("\n== unresolved groups (kept untrimmed — review/re-run) ==\n")
df |> filter(sd_ct > threshold) |> print(n = Inf)

# A single well, or two wells that disagree, can't be resolved -> set to NA so that bio rep drops out.
# (0_drop rep 2: RNA18S1 wells 15.8 / 35.2, DUSP11 21.9 / 29.7, single ISG15 well 34.7 -> failed sample)
df <- df |> mutate(unresolvable = n_wells == 1 | (n_wells == 2 & sd_ct > threshold))

cat("\n== unresolvable groups set to NA ==\n")
df |> filter(unresolvable) |> print(n = Inf)

df <- df |>
  mutate(mean_ct = if_else(unresolvable, NA_real_, mean_ct)) |>
  select(-unresolvable) |>
  select(-sd_ct, -n_wells)

df_wide <- df |>
  pivot_wider(names_from = Target.Name, values_from = mean_ct) |>
  separate_wider_delim(Sample.Name, "_", names = c("dose", "delivery"), cols_remove = FALSE) |>
  mutate(dose = factor(dose, levels = c("0", "0.1", "0.5", "1")),
         delivery = factor(delivery, levels = c("transfect", "drop"))) |>
  # A sample (bio rep) missing any gene is a failed sample -> drop it from every analysis,
  # so all reference choices use the same set of samples
  filter(if_all(c(RNA18S1, PPIA, DUSP11, ISG15), \(x) !is.na(x)))

# Sanity check: the two reference genes should move together across samples
cat("\n== reference genes per sample (mean CT) ==\n")
df_wide |>
  summarize(across(c(RNA18S1, PPIA), mean), .by = c(delivery, dose)) |>
  mutate(diff = PPIA - RNA18S1) |>
  arrange(delivery, dose) |>
  print()


# Full ddCT -> stats -> plot pipeline for a given reference CT
# ref_ct: function(data) returning the reference CT per row
run_analysis <- function(ref_ct, ref_label, outfile) {
  cat("\n\n############ Reference:", ref_label, "############\n")

  df_dct <- df_wide |>
    mutate(ref_ct = ref_ct(pick(everything())),
           dct_DUSP11 = DUSP11 - ref_ct,
           dct_ISG15 = ISG15 - ref_ct)

  # ddCT: each delivery normalized to its own 0 dose (mean dCT of the 0 bio reps)
  control_avg <- df_dct |>
    filter(dose == "0") |>
    summarize(across(c(dct_DUSP11, dct_ISG15), \(x) mean(x, na.rm = TRUE), .names = "ctrl_{.col}"), .by = delivery)

  df_ddct <- df_dct |>
    left_join(control_avg, by = "delivery") |>
    mutate(ddct_DUSP11 = dct_DUSP11 - ctrl_dct_DUSP11,
           ddct_ISG15 = dct_ISG15 - ctrl_dct_ISG15,
           foldchange_DUSP11 = 2^(-ddct_DUSP11),
           foldchange_ISG15 = 2^(-ddct_ISG15)) |>
    select(-starts_with("ctrl_")) |>
    arrange(delivery, dose, Task)

  # Long format: one row per sample x replicate x gene
  df_long <- df_ddct |>
    select(delivery, dose, Task, ddct_DUSP11, ddct_ISG15) |>
    pivot_longer(starts_with("ddct_"), names_to = "gene", names_prefix = "ddct_", values_to = "ddct") |>
    filter(!is.na(ddct))

  # Stats on ddCT (log2 scale), per gene x delivery: ANOVA, assumption checks, Dunnett vs 0 dose
  cat("\n== ANOVA ==\n")
  df_long |> group_by(gene, delivery) |> anova_test(ddct ~ dose) |> print()

  cat("\n== Shapiro on ANOVA residuals ==\n")
  df_long |>
    group_by(gene, delivery) |>
    group_modify(~ broom::tidy(shapiro.test(residuals(aov(ddct ~ dose, data = .x)))) |> select(statistic, p.value))
  cat("\n== Levene ==\n")
  df_long |> group_by(gene, delivery) |> levene_test(ddct ~ dose) |> print()

  set.seed(42)
  stat_test <- df_long |>
    group_by(gene, delivery) |>
    dunnett_test(ddct ~ dose, ref.group = "0") |>
    add_significance("p.adj")
  cat("\n== Dunnett vs 0 dose ==\n")

  # Plot data: mean and SD on ddCT, back-transformed to fold change (asymmetric bounds)
  df_plot <- df_long |>
    summarize(mean_ddct = mean(ddct), sd_ddct = sd(ddct), n = n(), .by = c(gene, delivery, dose)) |>
    mutate(mean_fc = 2^(-mean_ddct),
           fc_lower = 2^(-(mean_ddct + sd_ddct)),
           fc_upper = 2^(-(mean_ddct - sd_ddct)))

  # Bar chart: fold change vs 0 dose, stars from Dunnett on ddCT
  stat_plot <- stat_test |>
    left_join(df_plot |> select(gene, delivery, group2 = dose, fc_upper), by = c("gene", "delivery", "group2")) |>
    mutate(y.position = fc_upper * 1.05)

  p <- ggplot(df_plot, aes(x = dose, y = mean_fc)) +
    geom_col(fill = "steelblue") +
    geom_errorbar(aes(ymin = fc_lower, ymax = fc_upper), width = 0.2) +
    ggpubr::stat_pvalue_manual(stat_plot, x = "group2", label = "p.adj.signif", hide.ns = FALSE) +
    facet_wrap(gene ~ delivery, scales = "free_y", nrow = 2, labeller = label_wrap_gen(multi_line = FALSE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "DUSP11 and ISG15 expression",
         subtitle = paste0("Normalized to ", ref_label, "; fold change vs 0 dose within each delivery method"),
         x = "Dose", y = "Fold change vs 0 dose") +
    theme_minimal()
  ggsave(outfile, plot = p, width = 8, height = 6, dpi = 300)

  invisible(df_plot)
}

# Reference = geometric mean of RNA18S1 and PPIA expression.
# Geo mean of 2^-CT values == 2^-(arithmetic mean of CTs), so average the CTs.
run_analysis(\(d) (d$RNA18S1 + d$PPIA) / 2, "geometric mean of RNA18S1 and PPIA",
             "TN072/72.5/TN072.5_foldchange.png")

# Reference = RNA18S1 only (as in 72.4), for comparison
run_analysis(\(d) d$RNA18S1, "RNA18S1 only",
             "TN072/72.5/TN072.5_foldchange_18S.png")

# Reference = PPIA only, for comparison
run_analysis(\(d) d$PPIA, "PPIA only",
             "TN072/72.5/TN072.5_foldchange_PPIA.png")
