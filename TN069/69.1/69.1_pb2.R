library(tidyverse)
library(dplyr)

df <- read.csv("TN069/69.1/TN069.1_PB2_filtered.csv")
df <- df %>%
  group_by(Sample.Name, Target.Name) %>%
  summarize(CT = mean(CT, na.rm = FALSE), .groups = "drop") %>%
  separate(Sample.Name, into = c("base_name", "timepoint", "type"), sep = "_") %>%
  pivot_wider(names_from = type, values_from = CT) %>%
  # calculate fold change using input dilution factor
  mutate(fc = 0.2 * 2^(input - enrich)) %>%
  mutate(log2fc = log2(fc)) %>%
  mutate(percent_input = fc * 100)
write.csv(df, "TN069/69.1/69.1_PB2_results.csv", row.names = FALSE)

# plot results
p <- ggplot(df, aes(x = timepoint, y = percent_input, color = base_name)) +
  geom_point() +
  theme_classic() +
  labs(title = "PB2 Recovery",
       x = "Incubation Duration (hours)",
       y = "Percent Input"
  )
print(p)
ggsave("TN069/69.1/69.1_PB2_percent_input.png", plot = p, width = 8, height = 5, dpi = 300)
