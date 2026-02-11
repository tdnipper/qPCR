library(tidyverse)

df <- read.csv("TN069/69.4/TN069.4_trizol_filtered.csv")

data <- df %>%
  group_by(Sample.Name, Target.Name, Task) %>%
  summarize(
    mean.CT = mean(CT, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  separate(Sample.Name, into = c("Sample", "enrichment"), sep = "_")
data <- data %>%
  pivot_wider(names_from = Sample, values_from = mean.CT) %>%
  mutate(dCT = trizol - input + log2(0.1^-1)) %>%
  mutate(foldchange = 2^(-dCT)) %>%
  mutate(percent = 100 * foldchange)
write.csv(data, "TN069/69.4/TN069.4_trizol_percent_input.csv", row.names = FALSE)
p <- ggplot(data, aes(x = enrichment, y = percent)) +
    geom_point(position = position_dodge(width = 0.5)) +
    labs(
      title = "Trizol Elution",
      y = "Percent Input"
    ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
print(p)
# ggsave("TN069/69.4/TN069.4_trizol_percent_input_plot.png", plot = p, width = 8, height = 6, dpi = 300)