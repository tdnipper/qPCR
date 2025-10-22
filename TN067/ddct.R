library(readxl)
library(dplyr)
library(tidyr)
library(writexl)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(car)

df <- read_excel("TN067/TN067.1_filtered.xlsx")

# Calculate mean CT, pivot wider, and convert to numeric
data <- df %>%
    group_by(`Sample Name`, `Target Name`, `Task`) %>%
    summarise(CT = mean(CT, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = `Target Name`, values_from = CT)

target_cols <- setdiff(names(data), c("Sample Name", "Task"))
data[target_cols] <- lapply(data[target_cols], as.numeric)

# Calculate delta Ct for all targets
data <- data %>%
    group_by(`Sample Name`, `Task`) %>%
    mutate(across(
        .cols = all_of(target_cols),
        .fns = ~ .x - RNA18S1,
        .names = "delta_ct_{.col}"
    )) %>%
    ungroup()

# Pivot to long format for delta Ct
delta_long <- data %>%
    select(`Sample Name`, `Task`, starts_with("delta_ct_")) %>%
    pivot_longer(
        cols = starts_with("delta_ct_"),
        names_to = "Target",
        values_to = "delta_ct",
        names_prefix = "delta_ct_"
    )

# Get reference delta Ct for Control only
ref_delta <- delta_long %>%
    filter(`Sample Name` == "mock") %>%
    select(Target, delta_ct, Task) %>%
    rename(ref_delta_ct = delta_ct)

# Join reference and calculate ddCT and fold change
fold_change_data <- delta_long %>%
    left_join(ref_delta, by = c("Target", "Task")) %>%
    mutate(ddct = delta_ct - ref_delta_ct,
           fold_change = 2^(-ddct)) %>%
    filter(!`Sample Name` %in% c("mock", NA), `Target` != "RNA18S1")

# Export to Excel
write_xlsx(fold_change_data, "TN067/ddct_results.xlsx")

# canonicalize sample name (remove space), drop NA samples and remove 'mock'
fold_change_data <- fold_change_data %>%
  rename(Sample_Name = `Sample Name`) %>%
  filter(!is.na(Sample_Name), !Sample_Name %in% c("mock"))

# --- NORMALITY TESTS (Shapiro–Wilk per group) ---
by(fold_change_data$ddct, fold_change_data$Sample_Name, shapiro.test)
normality <- fold_change_data %>%
  group_by(Sample_Name) %>%
  summarise(
    W = shapiro.test(ddct)$statistic,
    p_value = shapiro.test(ddct)$p.value
  )
print(normality)

# --- HOMOGENEITY OF VARIANCE (Levene’s test) ---
levene <- leveneTest(ddct ~ Sample_Name, data = fold_change_data)
print(levene)

# Compute and adjust p-values
pvals <- fold_change_data %>%
  pairwise_wilcox_test(ddct ~ Sample_Name, ref.group = "T0", p.adjust.method = "holm")

pvals <- pvals %>%
  arrange(group2) %>%
  mutate(y.position = max(fold_change_data$fold_change, na.rm = TRUE) * (1 + 0.05 * row_number()))

# Plot
ymin <- 0
ymax <- max(fold_change_data$fold_change, na.rm = TRUE)
xorder <- c("T0", "T8", "T24", "T48")
p <- ggplot(fold_change_data, aes(x = factor(Sample_Name, levels = xorder), y = fold_change)) +
    geom_jitter(position = position_jitterdodge()) +
    coord_cartesian(ylim = c(ymin, ymax * 1.2)) +
    labs(title = "DUSP11 mRNA during multicycle infection",
         x = "Hours post-infection",
         y = "Fold Change / Mock (18S)") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16)
  ) +
  stat_pvalue_manual(pvals, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE)

# Save the plot
ggsave("TN067/fold_change_plot.png", plot = p, width = 8, height = 6)