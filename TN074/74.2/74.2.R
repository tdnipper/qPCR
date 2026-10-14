library(tidyverse)
library(rstatix)

df <- read_csv("TN074/74.2/TN074.2_DUSP11_filtered.csv", name_repair = "universal") 

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

# Calculate dCT
df_dct <- df |>
  pivot_wider(names_from = Target.Name, values_from = mean_ct) |>
  mutate(dct= DUSP11 - RNA18S1)

# Split into dox and no dox
df_dox <- df_dct |>
  filter(!str_detect(Sample.Name, pattern = "nodox")) |>
  select(Sample.Name, dct)
control_avg_dox <- df_dox |>
  filter(Sample.Name == "dox_mock") |>
  summarize(avg_dct = mean(dct)) |>
  pull(avg_dct)
df_dox_ddct <- df_dox |>
  mutate(ddct = dct - control_avg_dox)

df_nodox <- df_dct |>
  filter(str_detect(Sample.Name, pattern = "nodox")) |>
  select(Sample.Name, dct)
control_avg_nodox <- df_nodox |>
  filter(Sample.Name == "nodox_mock") |>
  summarize(avg_dct = mean(dct)) |>
  pull(avg_dct)
df_nodox_ddct <- df_nodox |>
  mutate(ddct = dct - control_avg_nodox)

# join dox and no dox
df_ddct <- bind_rows(df_dox_ddct, df_nodox_ddct)
df_ddct <- df_ddct |>
  mutate(foldchange = 2^(-ddct))

# Stats on ddCT during infection
# Welch t-tests on the log2 scale (ddCT), infected vs mock within each dox condition.
# dox and no dox are not cross comparable, so no ANOVA across all four groups.
stat_test <- df_ddct |>
  mutate(
    dox = if_else(str_detect(Sample.Name, "nodox"), "nodox", "dox"),
    infection = str_extract(Sample.Name, "mock|infected")
  ) |>
  group_by(dox) |>
  t_test(ddct ~ infection, ref.group = "mock", var.equal = FALSE) |>
  add_significance("p") |>
  mutate(group1 = paste(dox, group1, sep = "_"), group2 = paste(dox, group2, sep = "_"))
print(stat_test)

# p-values go into the exported csv (on the infected rows; mock rows are the reference)
df_ddct_out <- df_ddct |>
  left_join(
    stat_test |> select(Sample.Name = group2, t_statistic = statistic, df, p, p.signif),
    by = "Sample.Name"
  )
write_csv(df_ddct_out, "TN074/74.2/TN074.2_DUSP11_ddct.csv")
write_csv(stat_test, "TN074/74.2/TN074.2_DUSP11_ttests.csv")

df_plot <- df_ddct |>
  select(Sample.Name, ddct, foldchange) |>
  mutate(mean_fc = mean(foldchange), sd = sd(foldchange), .by = Sample.Name) |>
  select(Sample.Name, mean_fc, sd) |>
  distinct()
df_plot <- df_plot |> 
  mutate(Sample.Name = factor(Sample.Name, levels = c("dox_mock", "dox_infected", "nodox_mock", "nodox_infected")))

# Bracket positions sit above the taller bar+SD of each pair
stat_plot <- stat_test |>
  left_join(df_plot |> transmute(group1 = as.character(Sample.Name), top1 = mean_fc + sd), by = "group1") |>
  left_join(df_plot |> transmute(group2 = as.character(Sample.Name), top2 = mean_fc + sd), by = "group2") |>
  mutate(y.position = pmax(top1, top2) * 1.08)

p <- ggplot(df_plot, aes(x = Sample.Name, y = mean_fc)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_fc - sd, ymax = mean_fc + sd), width = 0.2) +
  ggpubr::stat_pvalue_manual(stat_plot, label = "p.signif", tip.length = 0.02, hide.ns = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "DUSP11 during infection ", x = "Sample Name", y = "Foldchange") +
  theme_minimal()
ggsave("TN074/74.2/TN074.2_DUSP11_foldchange.png", plot = p, width = 6, height = 4, dpi = 300)
