library(tidyverse)

data1 <- read.csv("TN066/66.10/TN066_10_enrichment_foldchange.csv")
data2 <- read.csv("TN066/66.11/TN066_11_enrichment_foldchange.csv")

data1 <- data1 %>%
  select(Target.Name, infection, FAB, V5)

data2 <- data2 %>%
  select(Target.Name, infection, FAB, V5)

data_combined <- data1 %>%
  inner_join(data2, by = c("Target.Name", "infection"), suffix = c("_10", "_11"))

data_combined <- data_combined %>%
  mutate(FAB = FAB_10 - FAB_11,
		 V5 = V5_10 - V5_11)

write_csv(data_combined %>% select(-FAB_10, -FAB_11, -V5_10, -V5_11), "TN066/66.11/66.10_11_combined_enrichment_foldchange.csv")