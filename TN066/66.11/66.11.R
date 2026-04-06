library(tidyverse)

data <- read.csv("TN066/66.11/TN066.11_filtered.csv") %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(mean_ct = mean(CT, na.rm = TRUE), .groups = 'drop')

data <- data %>%
  separate(Sample.Name, into = c("infection", "resin"), sep = "_")

data <- data %>%
  pivot_wider(names_from = resin, values_from = mean_ct) %>%
  mutate(enrich = FAB - V5) %>%
  mutate(enrich_foldchange = 2^-enrich)

write_csv(data, "TN066/66.11/TN066_11_enrichment_foldchange.csv")

p <- ggplot(data, aes(x = Target.Name, y = enrich_foldchange, color = infection)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  theme_classic() +
  labs(title = "RIG-I RIP", x = "", y = "Fold Change (FAB - V5)", color = "", caption = "66.11")
ggsave("TN066/66.11/TN066_11_enrichment_foldchange.png", plot = p, width = 10, height = 6, dpi=300)