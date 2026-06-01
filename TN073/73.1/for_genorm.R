library(tidyverse)

data <- read_csv("TN073/73.1/TN073.1_filtered.csv")

data <- data |>
  group_by(`Sample Name`, `Target Name`, `Task`) |>
  summarize(mean_CT = mean(CT), .groups = "drop")

data_genorm <- data |>
  pivot_wider(names_from = `Target Name`, values_from = mean_CT) |>
  select(-c(`Sample Name`, `Task`))

write_csv(data_genorm, "TN073/73.1/TN073.1_genorm.csv")