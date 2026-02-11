library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
library(ggpubr)

df <- read_excel("TN064/TN064.1_filtered.xlsx")

# Calculate mean CT, pivot wider, and convert to numeric
data <- df %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(CT = mean(CT, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = `Target Name`, values_from = CT)

data_long <- data %>%
  pivot_longer(cols = -`Sample Name`, names_to = "Target", values_to = "CT") %>%
  pivot_wider(names_from = `Sample Name`, values_from = CT) %>%
  mutate(across(.cols = -c("Target", "input"),
                .fns = ~ .x - input,
                .names = "delta_ct_{.col}")) %>%
  select(Target, starts_with("delta_ct_")) %>%
  mutate(across(.cols = starts_with("delta_ct_"),
                .fns = ~ 2^-.x,
                .names = "foldchange_{.col}")) %>%
  rename_with(~ sub("delta_ct_", "", .x), starts_with("foldchange_")) %>%
  select(Target, starts_with("foldchange_")) %>%
  pivot_longer(cols = -Target,
               names_to = "Capture",
               values_to = "Delta_CT",
               names_prefix = "foldchange_capture")


write_xlsx(data_long, "TN064/dct_results.xlsx")

# Graph dct data
p <- ggplot(data_long, aes(x = Target, y = Delta_CT, color = Capture)) +
  geom_point(size = 2, alpha = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Depletion Efficiency (Capture - Input)",
    x = "Target Gene",
    y = "foldchange",
    color = "Capture Oligo") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave("TN064/dct_plot.png", plot = p, width = 10, height = 6, dpi = 300)