library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
library(ggpubr)

df <- read_excel("TN066/TN066.1_filtered.xlsx")

# Calculate mean CT, pivot wider, and convert to numeric
data <- df %>%
    group_by(`Sample Name`, `Target Name`) %>%
    summarise(CT = mean(CT, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = `Sample Name`, values_from = CT)
# Calculate the percent recovery for each enriched sample
# correcting for dilution factor
data_dct <- data %>%
  mutate(across(.cols = `infected enriched`,
                .fns  = ~ . - `infected input` - log2(0.02^-1),
                .names = "delta_ct_infected")) %>%
  mutate(across(.cols = `mock enriched`,
                .fns  = ~ . - `mock input` - log2(0.04^-1),
                .names = "delta_ct_mock")) %>%
  select(`Target Name`, starts_with("delta_ct_")) %>%
  mutate(across(.cols = starts_with("delta_ct_"),
                .fns = ~ 2^-.x * 100,
                .names = "percent_recovery_{.col}")) %>%
  rename_with(~ sub("delta_ct_", "", .x), starts_with("percent_recovery_")) %>%
  select(`Target Name`, starts_with("percent_recovery_")) %>%
  pivot_longer(cols = -`Target Name`,
               names_to = "Condition",
               values_to = "percent_recovery",
               names_prefix = "percent_recovery_")

write_xlsx(data_dct, "TN066/rip_results.xlsx")

data_alt <- data %>%
  mutate(dct = `infected enriched` - `mock enriched`) 

rna18s1_dct <- data_alt %>%
  filter(`Target Name` == "RNA18S1") %>%
  pull(dct)

data_alt <- data_alt %>%
  mutate(dct_norm = dct - rna18s1_dct) %>%
  select(`Target Name`, dct, dct_norm)
write_xlsx(data_alt, "TN066/rip_dct_normalized.xlsx")

# Plot percent recovery data
p <- ggplot(data_dct, aes(x = `Target Name`, y = percent_recovery, color = Condition)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "NP Recovery Infected vs Mock",
    x = "Target Gene",
    y = "Percent Recovery",
    color = "Status") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave("TN066/percent_recovery_plot.png", plot = p, width = 10, height = 6, dpi = 300)

# Plot dct normalized data
p2 <- ggplot(data_alt, aes(x = `Target Name`, y = dct_norm)) +
  geom_point(size = 2, alpha = 0.7, color = "#253755") +
  labs(title = "NP Enrichment DCT Infected - Mock",
    x = "Target Gene",
    y = "DCT - 18S") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave("TN066/dct_normalized_plot.png", plot = p2, width = 10, height = 6, dpi = 300)

p3 <- ggplot(data_alt, aes(x = `Target Name`, y = dct)) +
  geom_point(size = 2, alpha = 0.7, color = "#253755") +
  labs(title = "NP Enrichment DCT Infected - Mock",
    x = "Target Gene",
    y = "DCT") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
ggsave("TN066/dct_plot.png", plot = p3, width = 10, height = 6, dpi = 300)