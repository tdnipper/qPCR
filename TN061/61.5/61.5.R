library(dplyr)
library(tidyverse)

df <- read.csv("TN061/61.5/TN061.5_filtered.csv") %>%
  group_by(Sample.Name, Target.Name) %>%
  summarize(CT = mean(CT, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Sample.Name, values_from = CT)

# Drop FFLUC row because it didn't amplify
df <- df %>% filter(Target.Name != "FFLUC")

# Calculate enrichment from input for each sample
df_percent <- df %>%
  mutate(`+RnaseH` = `+RnaseH_enrich` - `+RnaseH_input` - log2(0.1^-1)) %>%
  mutate(`+RnaseH_percent_input` = 2^(-`+RnaseH`)) %>%
  mutate(`-RnaseH` = `-RnaseH_enrich` - `-RnaseH_input` - log2(0.1^-1)) %>%
  mutate(`-RnaseH_percent_input` = 2^(-`-RnaseH`)) %>%
  mutate(foldchange = `+RnaseH_percent_input` / `-RnaseH_percent_input`) %>%
  select(Target.Name, `+RnaseH_percent_input`, `-RnaseH_percent_input`, foldchange)
write.csv(df_percent, "TN061/61.5/percent_input.csv", row.names = FALSE)

p <- ggplot(df_percent, aes(x = Target.Name, y = foldchange)) +
  geom_point() + 
  theme_classic() +
  labs(title = "Fold Change (+RnaseH / -RnaseH)",
       x = "Target Gene",
       y = "Fold Change")
print(p)
ggsave("TN061/61.5/foldchange_plot.png", plot = p, width = 8, height = 6, dpi = 300)