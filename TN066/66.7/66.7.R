library(tidyverse)

df <- read.csv("TN066/66.7/TN066.7_filtered.csv")

data <- df %>%
    group_by(Sample.Name, Target.Name, Task) %>%
    summarise(CT = mean(CT, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Sample.Name, values_from = CT)

data_enrich <- data %>%
    mutate(mock = mock_FLAG - mock_IgG,
        infect = infected_FLAG - infected_IgG) %>%
    mutate(across(.cols = c(mock, infect),
        ~ 2^-.x))

data_enrich_long <- data_enrich %>%
    select(Target.Name, Task, mock, infect) %>%
    pivot_longer(cols = c(infect, mock),
        names_to = "Condition",
        values_to = "normalized_foldchange",
        names_prefix = "foldchange_")

p <- ggplot(data_enrich_long, aes(x = Target.Name, y = normalized_foldchange, color = Condition)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    theme_classic()
ggsave("TN066/66.7/66.7_foldchange.png", plot = p, width = 10, height = 6, dpi=300)
