library(tidyverse)

data <- read.csv("TN066/66.10/TN066.10_filtered.csv") %>%
  group_by(Sample.Name, Target.Name) %>%
  summarize(mean_ct = mean(CT, na.rm = TRUE), .groups = 'drop')

data_eluate <- data %>%
  filter(!grepl("input", Sample.Name, ignore.case = TRUE)) %>%
  separate(Sample.Name, into = c("condition", "infection", "ab"), sep = "_") %>%
  pivot_wider(names_from = ab, values_from = mean_ct) %>%
  mutate(enrich = FAB - V5) %>%
  mutate(enrich_foldchange = 2^-enrich) 

write.csv(data_eluate, "TN066/66.10/TN066_10_enrichment_foldchange.csv")

p <- ggplot(data_eluate, aes(x = Target.Name, y = enrich_foldchange, color = infection)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  theme_classic() +
  labs(title = "RIG-I RIP", x = "", y = "Fold Change (FAB - V5)", color = "", caption = "66.10")
ggsave("TN066/66.10/TN066_10_enrichment_foldchange.png", plot = p, width = 10, height = 6, dpi=300)