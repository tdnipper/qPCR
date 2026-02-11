library(tidyverse)

df <- read.csv("TN069/69.3/TN069.3_biotin_filtered.csv")

data <- df %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean_Ct = mean(CT, na.rm = TRUE), .groups = "drop") %>%
  separate(Sample.Name, into = c("sample", "condition"), sep = "_") %>%
  pivot_wider(names_from = condition, values_from = mean_Ct) %>%
  mutate(across(.cols = c("flowthrough", "elute-heat", "elute+heat"),
                .fns = ~ .x - input + log2(0.1^-1),
                .names = "deltaCT_{.col}")) %>%
  mutate(across(.cols = c("deltaCT_flowthrough", "deltaCT_elute-heat", "deltaCT_elute+heat"),
                .fns = ~ 2^-.x,
                .names = "foldchange_{.col}")) %>%
  rename_with(~ sub("deltaCT_", "", .x), starts_with("foldchange_")) %>%
  mutate(across(.cols = c("foldchange_flowthrough", "foldchange_elute-heat", "foldchange_elute+heat"),
                .fns = ~ .x * 100,
                .names = "percent_{.col}")) %>%
  rename_with(~ sub("foldchange_", "", .x), starts_with("percent_"))

subset_data <- data %>%
  select(sample, Target.Name, Task,
         percent_flowthrough, `percent_elute-heat`, `percent_elute+heat`) %>%
  pivot_longer(cols = starts_with("percent_"),
               names_to = "condition",
               values_to = "percent_input",
               names_prefix = "percent_")

p <- ggplot(subset_data, aes(x = condition, y = percent_input, color = sample, group = sample)) +
  geom_point(position = position_dodge(width = 0.5)) +
  labs(title = "5'ppp FLUC Recovery",
       x = "Condition",
       y = "Percent input",
       color = "Biotin solvent",
       caption = "69.3, biotin solvent test") +
  theme_minimal()
print(p)
ggsave("TN069/69.3/TN069.3_biotin_percent_input_plot.png", plot = p, width = 8, height = 6, dpi = 300)