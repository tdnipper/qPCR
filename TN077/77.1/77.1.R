library(tidyverse)
library(rstatix)

df <- read_csv("TN077/77.1/TN077.1_filtered.csv", name_repair = "universal") 

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

# dCT vs 18S, then ddCT vs the mean 0uM dCT of each gene
ddct <- df |>
  pivot_wider(names_from = Target.Name, values_from = mean_ct) |>
  pivot_longer(-c(Sample.Name, Task, RNA18S1), names_to = "gene", values_to = "ct") |>
  mutate(dct = ct - RNA18S1,
         ddct = dct - mean(dct[Sample.Name == "0uM"]),
         Sample.Name = factor(Sample.Name, levels = c("0uM", "20uM", "50uM", "100uM")),
         .by = gene)

# Stats on ddCT, per gene: ANOVA, Levene, Dunnett of each Sample.Name vs 0uM
by_gene <- group_by(ddct, gene)
print(anova_test(by_gene, ddct ~ Sample.Name))
print(levene_test(by_gene, ddct ~ Sample.Name))

set.seed(24)
stat_test <- by_gene |>
  dunnett_test(ddct ~ Sample.Name, ref.group = "0uM") |>
  add_significance("p.adj")

# Paired bar chart: pre and total side by side per Sample.Name, one panel per amino acid.
# Mean and SD on ddCT, back-transformed to fold change vs 0uM (asymmetric bounds)
split_gene <- function(d) {
  d |> separate(gene, c("form", "aa"), sep = "_") |>
    mutate(form = factor(form, levels = c("pre", "total")))
}

stars <- stat_test |>
  filter(p.adj.signif != "ns") |>
  split_gene() |>
  select(aa, form, Sample.Name = group2, p.adj.signif) |>
  mutate(Sample.Name = factor(Sample.Name, levels = levels(ddct$Sample.Name)))

df_plot <- ddct |>
  split_gene() |>
  summarize(mean_ddct = mean(ddct), sd_ddct = sd(ddct), .by = c(aa, form, Sample.Name)) |>
  mutate(mean_fc = 2^(-mean_ddct),
         fc_lower = 2^(-(mean_ddct + sd_ddct)),
         fc_upper = 2^(-(mean_ddct - sd_ddct))) |>
  left_join(stars, by = c("aa", "form", "Sample.Name"))

dodge <- position_dodge(width = 0.8)
p <- ggplot(df_plot, aes(x = Sample.Name, y = mean_fc, fill = form)) +
  geom_col(position = dodge, width = 0.7) +
  geom_errorbar(aes(ymin = fc_lower, ymax = fc_upper), position = dodge, width = 0.25) +
  geom_text(aes(y = fc_upper, label = p.adj.signif), position = dodge, vjust = -0.4, na.rm = TRUE) +
  facet_wrap(~ aa) +
  scale_fill_manual(values = c(pre = "#2A3752", total = "#C4694F")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "", y = "Fold change", title = "tRNA (pre and total) change during pol III inhibition") +
  theme_minimal()
ggsave("TN077/77.1/TN077.1_foldchange.png", plot = p, width = 8, height = 4, dpi = 300)
