library(tidyverse)

df <- read.csv("TN069/69.2/TN069.2_filtered.csv")

data <- df %>%
  group_by(Sample.Name, Target.Name) %>%
  summarize(mean_Ct = mean(CT, na.rm = TRUE), .groups = "drop") %>%
  separate(Sample.Name, into = c("sample", "timepoint", "condition"), sep = "_") %>%
  pivot_wider(names_from = condition, values_from = mean_Ct) %>%
  mutate(delta_Ct = enrich - input + log2(0.1^-1)) %>%
  mutate(foldchange = 2^-(delta_Ct)) %>%
  mutate(percent = foldchange * 100)

p <- ggplot(data, aes(x = Target.Name, y = percent, color = timepoint, group = sample)) +
  geom_point(aes(shape = sample)) +
  labs(title = "Capping timecourse recovery",
       x = "Timepoint",
       y = "Percent input concentration") +
  theme_minimal()
print(p)
ggsave("TN069/69.2/TN069.2_percent_input_plot.png", plot = p, width = 8, height = 6, dpi = 300)