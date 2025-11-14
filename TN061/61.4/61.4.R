library(ggplot2)
library(readxl)
library(dplyr)
library(writexl)

df <- read.csv("TN061/TN061.4_filtered.csv")

# calculate mean CT for each Sample Name and Target Name
df <- df %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(CT = mean(CT, na.rm = TRUE), .groups = "drop")
write_xlsx(df, "TN061/61.4_ct_values.xlsx")

# Plot CT values
p <- ggplot(df, aes(x = `Target Name`, y = CT, color = `Sample Name`)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "RnaseH Treatment",
       x = "Sample Name",
       y = "CT",
       color = "Gene",
       caption = "No capping/enrichment") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave("TN061/61.4_ct_values_plot.png", plot = p, width = 10, height = 6, dpi = 300)
