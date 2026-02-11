library(tidyverse)

# Import filtered data
df <- read.csv("TN069/69.4/TN069.4_filtered.csv")

data <- df %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(
    mean.CT = mean(CT, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  separate(Sample.Name, into = c("Sample", "enrichment"), sep = "_") %>%
  pivot_wider(names_from = enrichment, values_from = mean.CT)

# calculate dCT, foldchange, and percent recovery
data <- data %>%
  mutate(across(
    .cols = c("elute", "flowthrough"),
    .fns = ~ .x - input + log2(0.1^-1),
    .names = "deltaCT_{col}"
  )) %>%
  mutate(across(
    .cols = starts_with("deltaCT_"),
    .fns = ~ 2 ^ (-.x),
    .names = "foldchange_{col}"
  )) %>%
  rename_with(
    .cols = starts_with("foldchange_"),
    .fn = ~ str_replace(.x, "foldchange_deltaCT_", "foldchange_")
  ) %>%
  mutate(across(
    .cols = starts_with("foldchange_"),
    .fns = ~ .x * 100,
    .names = "percent_input_{col}"
  )) %>%
  rename_with(
    .cols = starts_with("percent_input_"),
    .fn = ~ str_replace(.x, "percent_input_foldchange_", "percent_input_")
  )

# subset and pivot longer for plotting
subset_data <- data %>%
  select(Sample, Target.Name, Task, starts_with("percent_input_")) %>%
  rename_with(
    .cols = starts_with("percent_input_"),
    .fn = ~ str_replace(.x, "percent_input_", "")
  ) %>%
  pivot_longer(
    cols = c("flowthrough", "elute"),
    names_to = "enrichment",
    values_to = "percent_input"
  )

# import trizol data from other analysis
trizol_data <- read.csv("TN069/69.4/TN069.4_trizol_percent_input.csv") %>%
  select(Sample = enrichment, Target.Name, Task, percent) %>%
  mutate(enrichment = "trizol") %>%
  rename(percent_input = percent)

# join trizol data to subset_data
subset_data <- bind_rows(subset_data, trizol_data)

p <- ggplot(subset_data, aes(x = enrichment, y = percent_input, color = Sample)) +
    geom_point(position = position_dodge(width = 0.5)) +
    labs(
      title = "Percent Input +- Rnasin",
      y = "Percent Input"
    ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
print(p)
ggsave("TN069/69.4/TN069.4_percent_input_plot.png", plot = p, width = 8, height = 6, dpi = 300)