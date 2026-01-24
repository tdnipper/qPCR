library(tidyverse)

df <- read.csv("TN069/69.5/TN069.5_filtered.csv")

data <- df %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(
    Mean.CT = mean(CT, na.rm = TRUE), 
    .groups = 'drop'
    ) %>%
  separate(Sample.Name, into = c("target", "enrich"), sep = "_") %>%
  pivot_wider(names_from = enrich, values_from = Mean.CT)
  
# calculate ddCT, foldchange, and percent recovery
data <- data %>%
  mutate(
    across(
      .cols = c("enrich", "flowthrough"),
      .fns = ~ .x - input  + log2(0.1^-1),
      .names = "dCT_{.col}" 
    )
  ) %>%
  mutate(across(
    .cols = starts_with("dCT_"),
    .fns = ~ 2^(-.x),
    .names = "foldchange_{.col}"
  )) %>%
  rename_with(
    .cols = starts_with("foldchange_"),
    .fn = ~ str_replace(.x, "foldchange_dCT_", "foldchange_")
  ) %>%
  mutate(
    across(
      .cols = starts_with("foldchange_"),
      .fns = ~ .x * 100,
      .names = "percent_input_{.col}"
    )) %>%
  rename_with(
    .cols = starts_with("percent_input_"),
    .fn = ~ str_replace(.x, "percent_input_foldchange_", "percent_input_")
  )
write.csv(data, "TN069/69.5/TN069.5_recovery_calculations.csv", row.names = FALSE)

# plot percent recovery
plot_data <- data %>%
  select(Target.Name, starts_with("percent_input_"), target) %>%
  rename_with(
    .cols = starts_with("percent_input_"),
    .fn = ~ str_replace(.x, "percent_input_", "")
  ) %>%
  pivot_longer(
    cols = c("enrich", "flowthrough"),
    names_to = "condition",
    values_to = "percent_input"
  )

p <- ggplot(plot_data, aes(x = condition, y = percent_input, color = Target.Name, shape = target)) +
  geom_point(stat = "identity", position = position_dodge(width = 0.5)) +
  labs(
    title = "Percent Input Recovery for TN069.5",
    x = "Target Name",
    y = "Percent Input Recovery (%)",
    color = "",
    shape = ""
  ) +
  theme_minimal()
print(p)
ggsave("TN069/69.5/TN069.5_percent_input.png", plot = p, width = 8, height = 6, dpi = 300)