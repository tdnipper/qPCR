library(tidyverse)

data <- read_tsv("TN073/73.1/genorm_results.tsv")

p <- ggplot(data, aes(x = reorder(Gene, -Rank), y = `M-value`)) +
  geom_line(group = 1, color = "#2A3752") +
  geom_point(color = "#2A3752") +
  labs(title = "Gene Stability Analysis (geNorm)",
	   x = "",
	   y = "M-value") +
  theme_classic() +
  theme(legend.title = element_blank())

ggsave("TN073/73.1/genorm_stability_plot.png", plot = p, width = 8, height = 6, dpi=300)