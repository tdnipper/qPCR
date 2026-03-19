library(tidyverse)

data <- read.csv("TN066/66.9/TN066.9_input.csv") %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean_ct = mean(CT, na.rm = TRUE), .groups = 'drop') %>%
  mutate(mean_ct = mean_ct - log2(.044^-1))

# Create filter to grab original data
samples = data %>%
  select(Sample.Name) %>%
  distinct() %>%
  pull()

data_enrich_mock <- read.csv("TN066/66.9/TN066_9_mock.csv") %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean_ct = mean(CT, na.rm = TRUE), .groups = 'drop') %>%
  filter(Sample.Name %in% samples)

data_enrich_infect <- read.csv("TN066/66.9/TN066_9_infect.csv") %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean_ct = mean(CT, na.rm = TRUE), .groups = 'drop') %>%
  filter(Sample.Name %in% samples)

combined_data <- bind_rows(data_enrich_mock, data_enrich_infect)
full_data <- left_join(data, combined_data, by = c("Sample.Name", "Target.Name", "Task"), suffix = c("_input", "_enrich")) %>%
  mutate(norm_input = mean_ct_input - mean_ct_enrich) %>%
  select(Sample.Name, Target.Name, Task, norm_input)

p <- ggplot(full_data, aes(x = Target.Name, y = norm_input, color = Sample.Name, shape = Task)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  theme_classic() +
  labs(title = "Normalized CT (Input - Enrichment)", x = "Target", y = "Normalized CT")
ggsave("TN066/66.9/TN066_9_input_normalized.png", plot = p, width = 10, height = 6, dpi=300)