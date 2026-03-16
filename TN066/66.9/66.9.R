library(tidyverse)

# Join mock and infected data
mock <- read.csv("TN066/66.9/TN066_9_mock.csv")
infect <- read.csv("TN066/66.9/TN066_9_infect.csv")

data <- rbind(mock, infect)

# Calculate mean and standard deviation for each group
summary_data <- data %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean_ct = mean(CT, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Sample.Name, values_from = mean_ct)

data_enrich <- summary_data %>%
  mutate(mock_M2_norm = mock_M2 - mock_IgG,
		 infect_M2_norm = infect_M2 - infect_IgG) %>%
  mutate(across(.cols = c(mock_M2_norm, infect_M2_norm),
				~ 2^-.x))

data_enrich_long <- data_enrich %>%
  select(Target.Name, Task, mock_M2_norm, infect_M2_norm) %>%
  pivot_longer(cols = c(infect_M2_norm, mock_M2_norm),
			   names_to = "Condition",
			   values_to = "normalized_foldchange",
			   names_prefix = "foldchange_") %>%
  filter(Task != "-Dox")
write.csv(data_enrich_long, "TN066/66.9/TN066_9_combined_foldchange.csv", row.names = FALSE)

p <- ggplot(data_enrich_long, aes(x = Target.Name, y = normalized_foldchange, color = Condition)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  theme_classic() +
  ylim(0, 50)
ggsave("TN066/66.9/TN066_9_combined_foldchange.png", plot = p, width = 10, height = 6, dpi=300)

data_induction <- data %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean = mean(CT, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Task, values_from = mean) %>%
  mutate(norm_1 = `1` - `-Dox`) %>%
  mutate(norm_2 = `2` - `-Dox`) %>%
  pivot_longer(cols = c(norm_1, norm_2), names_to = "Condition", values_to = "normalized_CT")

p2 <- ggplot(data_induction, aes(x = Target.Name, y = normalized_CT, color = Condition)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  theme_classic() +
  ylim(-5, 10) +
  labs(title = "CT Target - (-Dox) control", x = "Target", y = "Normalized CT")
ggsave("TN066/66.9/TN066_9_induction_normalized.png", plot = p2, width = 10, height = 6, dpi=300)