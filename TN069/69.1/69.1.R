library(dplyr)
library(tidyverse)
library(ggplot2)

df <- read.csv("TN069/69.1/TN069.1_filtered.csv")
df <- df %>%
  group_by(Sample.Name, Target.Name) %>%
  summarize(CT = mean(CT, na.rm = FALSE), .groups = "drop") %>%
  separate(Sample.Name, into = c("base_name", "type"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = type, values_from = CT) %>%
  mutate(percent_input = 100 * 0.2 * 2^(input - enrich)) %>% 
  mutate(fc = 0.2 * 2^(input - enrich)) %>%
  mutate(log2fc = log2(fc))
write.csv(df, "TN069/69.1/69.1_results.csv", row.names = FALSE)

# filter for FLUC only
fluc <- df %>%
  filter(Target.Name == "FFLUC") %>%
  select(base_name, percent_input)

# Plot FLUC percent input
p <- ggplot(fluc, aes(x = base_name, y = percent_input)) +
  geom_point() +
  theme_classic() +
  labs(title = "FLUC Recovery",
       x = "sample_h.p.i.",
       y = "Percent Input"
  )
print(p)
ggsave("TN069/69.1/69.1_fluc_percent_input.png", plot = p, width = 8, height = 5, dpi = 300)


