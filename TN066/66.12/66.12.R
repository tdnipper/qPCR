library(tidyverse)

# 1. Load the data
data <- read_csv("TN066/66.12/TN066.12_filtered.csv")

# 2. Take mean of tech reps
data <- data %>%
  separate(`Sample Name`, into = c("infection", "resin"), sep = "_") %>%
  group_by(infection, resin, `Target Name`, Task) %>%
  summarize(mean_ct = mean(CT, na.rm = TRUE), .groups = 'drop')

# 3. Pivot and calculate foldchange between FAB and GFP
data_enrich <- data %>%
  pivot_wider(names_from = resin, values_from = mean_ct) %>%
  mutate(enrich = FAB - GFP) %>%
  mutate(enrich_foldchange = 2^-enrich)

write_csv(data_enrich, "TN066/66.12/TN066_12_enrichment_foldchange.csv")

# 4. Plot
p <- ggplot(data_enrich, aes(x = `Target Name`, y = enrich_foldchange, color = infection)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  theme_classic() +
  labs(title = "RIG-I RIP", x = "", y = "Fold Change (FAB - GFP)", color = "", caption = "66.12")
ggsave("TN066/66.12/TN066_12_enrichment_foldchange.png", plot = p, width = 10, height = 6, dpi=300)