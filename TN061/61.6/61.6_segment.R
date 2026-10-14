library(tidyverse)

df <- read_csv("TN061/61.6/TN061.6_segment_filtered.csv", name_repair = "universal")

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

df_dct <- df |>
  pivot_wider(names_from = Sample.Name, values_from = mean_ct) |>
  mutate(dct = enriched - input + log2(0.1^-1)) |>
  mutate(fold_change = 2^(-dct)) |>
  select(Target.Name, Task, dct, fold_change)

data_plot <- df_dct |>
  select(Target.Name, fold_change) |>
  mutate(mean_fold_change = mean(fold_change),
         sd_fold_change = sd(fold_change),
         .by = Target.Name) |>
  select(Target.Name, mean_fold_change, sd_fold_change) |>
  distinct()

p <- ggplot(data = data_plot, aes(x = Target.Name, y = mean_fold_change)) +
  geom_col(fill = "#2A3752", alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_fold_change - sd_fold_change,
                    ymax = mean_fold_change + sd_fold_change),
                width = 0.2,
                color = "#2A3752") +
  labs(title = "vRNA recovery after RnaseH", x = "", y = "Average Recovery (enriched/input)", caption = "61.6_segment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
ggsave("TN061/61.6/TN061.6_vRNA_recov ery.png", p, width = 6, height = 4, dpi = 300)