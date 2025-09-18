library(readxl)
df <- read_excel("TN052.13/filtered_results.xlsx")

# Calculate ΔCt for each sample
library(dplyr)
ddct_data <- df %>%
  group_by(`Sample Name`, `Target Name`, `Task`) %>%
  summarise(CT = mean(CT, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = `Target Name`, values_from = CT) %>%
  mutate(delta_ct = DUSP11 - RNA18S1)

# Reference sample: first sample
ref_delta_ct <- mean(ddct_data$delta_ct[1:3], na.rm = TRUE)

# Calculate ddCT and fold change
ddct_data <- ddct_data %>%
  mutate(ddct = delta_ct - ref_delta_ct,
         fold_change = 2^(-ddct))

# Plot foldchange
library(ggplot2)
ggplot(ddct_data, aes(x = `Sample Name`, y = fold_change)) +
  geom_boxplot() +
  geom_point(size = 2) +
  labs(title = "DUSP11 mRNA During Early Infection",
       x = "Hours post-infection",
       y = "Fold Change") +
  theme_classic() +
  theme(axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

# Save the plot
ggsave("TN052.13/fold_change_plot.png", width = 8, height = 6)