library(tidyverse)

df <- read.csv("TN069/69.2/TN069.2_flowthrough_filtered.csv")

data <- df %>%
  group_by(Sample.Name, Target.Name) %>%
  summarize(mean_Ct = mean(CT, na.rm = TRUE), .groups = "drop") %>%
  separate(Sample.Name, into = c("sample", "timepoint", "condition"), sep = "_") %>%
  pivot_wider(names_from = condition, values_from = mean_Ct) %>%
  mutate(across (.cols = c("enrich", "flowthrough"),
                .fns = ~ .x - input + log2(0.1^-1),
                .names = "deltaCT_{.col}")) %>%
  mutate(across (.cols = c("deltaCT_enrich", "deltaCT_flowthrough"),
                .fns = ~ 2^-.x,
                .names = "foldchange_{.col}")) %>%
  rename_with(~ sub("deltaCT_", "", .x), starts_with("foldchange_")) %>%
  mutate(across (.cols = c("foldchange_enrich", "foldchange_flowthrough"),
                .fns = ~ .x * 100,
                .names = "percent_{.col}")) %>%
  rename_with(~ sub("foldchange_", "", .x), starts_with("percent_"))

p_flowthrough <- ggplot(data, aes(x = timepoint, y = percent_flowthrough, color = timepoint, group = sample)) +
  geom_point(aes(shape = sample)) +
  labs(title = "Capping timecourse recovery (flowthrough)",
       x = "Timepoint",
       y = "Percent input concentration") +
  theme_minimal()
print(p_flowthrough)
ggsave("TN069/69.2/TN069.2_flowthrough_percent_input_plot.png", plot = p_flowthrough, width = 8, height = 6, dpi = 300)
p_enrich <- ggplot(data, aes(x = timepoint, y = percent_enrich, color = timepoint, group = sample)) +
  geom_point(aes(shape = sample)) +
  labs(title = "Capping timecourse recovery (enrichment)",
       x = "Timepoint",
       y = "Percent input concentration") +
  theme_minimal()
print(p_enrich)
ggsave("TN069/69.2/TN069.2_enrich_percent_input_plot.png", plot = p_enrich, width = 8, height = 6, dpi = 300)