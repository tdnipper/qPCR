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
  mutate(`+RnaseH_percent_input` = 2^(-`+RnaseH`)*100) %>%
  mutate(`-RnaseH` = `-RnaseH_enrich` - `-RnaseH_input`- log2(0.1^-1)) %>%
  mutate(`-RnaseH_percent_input` = 2^(-`-RnaseH`)*100) %>%
  select(Target.Name, `+RnaseH_percent_input`, `-RnaseH_percent_input`)

df_foldchange <- df %>%
  mutate(`foldchange` = 2^-(`+RnaseH_enrich` - `-RnaseH_enrich`)) %>%
  select(Target.Name, foldchange)