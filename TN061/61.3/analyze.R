library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
library(ggpubr)

df <- read_excel("TN061/61.3/filtered_results.xlsx")

# Calculate mean CT, pivot wider, and convert to numeric
data <- df %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(CT = mean(CT, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = `Target Name`, values_from = CT)

target_cols <- setdiff(names(data), "Sample Name")
data[target_cols] <- lapply(data[target_cols], as.numeric)

# Calculate delta Ct for all targets
data <- data %>%
  mutate(across(
    .cols = target_cols,
    .fns = ~ .x - RNA18S1,
    .names = "delta_ct_{.col}"
  ))

# Pivot to long format for delta Ct
delta_long <- data %>%
  select(`Sample Name`, starts_with("delta_ct_")) %>%
  pivot_longer(
    cols = starts_with("delta_ct_"),
    names_to = "Target",
    values_to = "delta_ct",
    names_prefix = "delta_ct_"
  )
print(delta_long)
# Get reference delta Ct for -RNaseH only
ref_delta <- delta_long %>%
  filter(`Sample Name` == "-RNaseH") %>%
  select(Target, delta_ct) %>%
  rename(ref_delta_ct = delta_ct)

# Join reference and calculate ddCT and fold change
fold_change_data <- delta_long %>%
  left_join(ref_delta, by = "Target") %>%
  mutate(ddct = delta_ct - ref_delta_ct,
         fold_change = 2^(-ddct)) %>%
  filter(`Sample Name` != "-RNaseH", `Target` != "RNA18S1")

# Export to Excel
write_xlsx(fold_change_data, "TN061/61.3/foldchange.xlsx")

# Plot
p <- ggplot(fold_change_data, aes(x = Target, y = fold_change)) +
  geom_point(size = 2) +
  labs(title = "Foldchange after RNaseH treatment",
       x = "",
       y = "Fold Change (18S)") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

ggsave("TN061/61.3/fold_change_plot.png", plot = p, width = 8, height = 6)