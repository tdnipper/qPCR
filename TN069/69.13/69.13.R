library(tidyverse)

# Import data and summarize mean ct
df <- read.csv("TN069/69.13/TN069.13_filtered.csv") %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean_ct = mean(CT), .groups = "drop")

# Calculate dCT between input and elutions
df_dct <- df %>%
  pivot_wider(names_from = Sample.Name, values_from = mean_ct) %>%
  mutate(across(.cols = c("flowthrough", "elute", "formamide"),
	.fns = ~ . - input + log2(0.1^-1),
	.names = "dCT_{col}")) %>%
  mutate(across(.cols = starts_with("dCT_"), .fns = ~ 2^(-.x), .names = "foldchange_{col}")) %>%
  rename_with(.cols = starts_with("foldchange_"), .fn = ~ str_replace(.x, "foldchange_dCT_", "foldchange_")) %>%
  mutate(across(.cols = starts_with("foldchange_"), .fns = ~ .x * 100, .names = "percent_input_{col}")) %>%
  rename_with(.cols = starts_with("percent_input_foldchange_"), .fn = ~ str_replace(.x, "percent_input_foldchange_", "percent_input_"))

plot_data <- df_dct %>%
  select(Target.Name, starts_with("percent_input"), Task) %>%
  pivot_longer(cols = starts_with("percent_input_"), names_to = "condition", values_to = "percent_input") %>%
  mutate(condition = str_replace(condition, "percent_input_", ""))
write.csv(plot_data, "TN069/69.13/TN069.13_percent_input.csv", row.names = FALSE)

missing_percentage <- plot_data %>%
  group_by(Target.Name, Task) %>%
  summarize(missing_percent = 100 - sum(percent_input), .groups = "drop")

# Compute label y-positions: stagger above max point per target
max_y <- plot_data %>%
  group_by(Target.Name) %>%
  summarize(max_y = max(percent_input), .groups = "drop")

missing_percentage <- missing_percentage %>%
  left_join(max_y, by = "Target.Name") %>%
  group_by(Target.Name) %>%
  mutate(y_pos = max_y + 5 + (row_number() - 1) * 4) %>%
  ungroup()

p <- ggplot(plot_data, aes(x = Target.Name, y = percent_input, color = condition)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.5)) +
  geom_text(data = missing_percentage,
            aes(x = Target.Name, y = y_pos,
                label = paste0("Task ", Task, ": ", round(missing_percent, 1), "% missing")),
            color = "black", size = 3, inherit.aes = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "Percent Input Recovery TN069.13") +
  theme_minimal()

ggsave("TN069/69.13/TN069.13_percent_input_plot.png", plot = p, width = 8, height = 6)