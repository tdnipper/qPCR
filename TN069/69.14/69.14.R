library(tidyverse)

data <- read.csv("TN069/69.14/TN069.14_filtered.csv") %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean_ct = mean(CT), .groups = "drop") %>%
  separate(Sample.Name, into = c("Sample", "eluate"), sep = "_")

data_dct <- data %>%
  pivot_wider(names_from = eluate, values_from = mean_ct) %>%
  mutate(across(.cols = c("flowthrough", "enrich"), .fns = ~ . - input + log2(0.1^-1), .names = "dCT_{col}")) %>%
  mutate(across(.cols = starts_with("dCT_"), .fns = ~ 2^(-.x), .names = "foldchange_{col}")) %>%
  rename_with(.cols = starts_with("foldchange_"), .fn = ~ str_replace(.x, "foldchange_dCT_", "foldchange_")) %>%
  mutate(across(.cols = starts_with("foldchange_"), .fns = ~ .x * 100, .names = "percent_input_{col}")) %>%
  rename_with(.cols = starts_with("percent_input_foldchange_"), .fn = ~ str_replace(.x, "percent_input_foldchange_", "percent_input_"))

plot_data <- data_dct %>%
  select(Sample, Target.Name, starts_with("percent_input"), Task) %>%
  pivot_longer(cols = starts_with("percent_input_"), names_to = "condition", values_to = "percent_input") %>%
  mutate(condition = str_replace(condition, "percent_input_", ""))
write.csv(plot_data, "TN069/69.14/TN069.14_percent_input.csv", row.names = FALSE)

p <- ggplot(plot_data, aes(x = Target.Name, y = percent_input, color = condition, shape=Sample)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.5)) +
  labs(title = "Percent Input Recovery TN069.14") +
  ylim(0,100) +
  theme_minimal()
ggsave("TN069/69.14/TN069.14_percent_input_plot.png", plot = p, width = 8, height = 6)