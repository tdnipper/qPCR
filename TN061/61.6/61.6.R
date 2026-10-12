library(tidyverse)

df <- read.csv("TN061/61.6/TN061.6_filtered.csv") |>
  group_by(Sample.Name, Target.Name, Task) |>
  summarize(CT = mean(CT, na.rm = TRUE), .groups = "drop")

# split dfs into +- RnaseH
# filter out 5'/3' fragments because these are not checked for enrichment
df_plus <- df |>
  filter(str_detect(Sample.Name, "RnaseH"), !str_detect(Target.Name, "5'|3'")) |>
  pivot_wider(names_from = Sample.Name, values_from = CT)
df_minus <- df |>
  filter(!str_detect(Sample.Name, "RnaseH"), !str_detect(Target.Name, "5'|3'")) |>
  pivot_wider(names_from = Sample.Name, values_from = CT)

# Calculate percent recovery
# Input is 10% frac so adjust by log2(0.1^-1) 
df_plus <- df_plus |>
  mutate(percent_recovery = `enriched_RnaseH` - `input_RnaseH` - log2(0.1^-1)) |>
  mutate(percent_recovery = 2^(-percent_recovery)) |>
  group_by(Target.Name) |>
  mutate(avg_percent_recovery = mean(percent_recovery, na.rm = TRUE)) |>
  mutate(sd_percent_recovery = sd(percent_recovery, na.rm = TRUE)) |>
  select(Target.Name, avg_percent_recovery, sd_percent_recovery) |>
  distinct()
df_minus <- df_minus |>
  mutate(percent_recovery = `enriched_untreated` - `input_untreated` - log2(0.1^-1)) |>
  mutate(percent_recovery = 2^(-percent_recovery)) |>
  group_by(Target.Name) |>
  mutate(avg_percent_recovery = mean(percent_recovery, na.rm = TRUE)) |>
  mutate(sd_percent_recovery = sd(percent_recovery, na.rm = TRUE)) |>
  select(Target.Name, avg_percent_recovery, sd_percent_recovery) |>
  distinct()

# Combine dfs and plot
df_combined <- df_plus |>
  select(Target.Name, avg_percent_recovery, sd_percent_recovery) |>
  rename(avg_percent_recovery_plus = avg_percent_recovery, sd_percent_recovery_plus = sd_percent_recovery) |>
  left_join(df_minus |>
              select(Target.Name, avg_percent_recovery, sd_percent_recovery) |>
              rename(avg_percent_recovery_minus = avg_percent_recovery, sd_percent_recovery_minus = sd_percent_recovery),
            by = "Target.Name") |>
  mutate(ratio = avg_percent_recovery_plus / avg_percent_recovery_minus)

write.csv(df_combined, "TN061/61.6/61.6_percent_recovery.csv", row.names = FALSE)

p <- df_combined |>
  pivot_longer(c(avg_percent_recovery_plus, avg_percent_recovery_minus,
                 sd_percent_recovery_plus, sd_percent_recovery_minus),
               names_to = c("statistic", "treatment"),
               names_pattern = "(avg|sd)_percent_recovery_(plus|minus)") |>
  pivot_wider(names_from = statistic, values_from = value) |>
  mutate(treatment = recode(treatment,
                            plus = "RnaseH",
                            minus = "Untreated")) |>
  ggplot(aes(x = Target.Name, y = avg, fill = treatment)) +
  geom_col(position = "identity", alpha = 0.7) +
  geom_errorbar(aes(ymin = avg - sd, ymax = avg + sd),
                position = "identity", width = 0.2) +
  labs(x = "", y = "Average Percent Recovery", fill = "Treatment") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
ggsave("TN061/61.6/61.6_percent_recovery_plot.png", p, dpi = 300)

# check ratio of 5'/3' fragments
df_5_3 <- df |>
  filter(str_detect(Sample.Name, "input")) |>
  filter(str_detect(Target.Name, "5'|3'") | str_detect(Target.Name, "18S"))

ratio_5_3 <- df_5_3 |>
  pivot_wider(names_from = Sample.Name, values_from = CT) |>
  mutate(delta_delta_ct =
           (`input_RnaseH` - `input_RnaseH`[Target.Name == "RNA18S1"]) -
           (`input_untreated` - `input_untreated`[Target.Name == "RNA18S1"])) |>
  mutate(ratio = 2^(-delta_delta_ct)) |>
  group_by(Target.Name) |>
  mutate(avg_ratio = mean(ratio, na.rm = TRUE)) |>
  mutate(sd_ratio = sd(ratio, na.rm = TRUE))

p_ratio <- ratio_5_3 |>
  select(Target.Name, avg_ratio, sd_ratio) |>
  filter(Target.Name != "RNA18S1") |>
  distinct() |>
  ggplot(aes(x = Target.Name, y = avg_ratio)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_errorbar(aes(ymin = avg_ratio - sd_ratio, ymax = avg_ratio + sd_ratio),
                width = 0.2) +
  labs(x = "", y = "Average Ratio of 5'/3' Fragments") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
ggsave("TN061/61.6/61.6_5_3_ratio_plot.png", p_ratio, dpi = 300)
