library(tidyverse)
library(rstatix)
library(ggpubr)

data <- read_csv("TN072/72.1/TN072.1_filtered.csv") %>%
  rename(sample_name = `Sample Name`, target_name = `Target Name`)

data <- data %>%
  group_by(sample_name, target_name, Task) %>%
  summarize(mean_ct = mean(CT), .groups = "drop")

# Calculate dCT of DUSP11 - 18S for each sample
data_dct <- data %>%
  pivot_wider(names_from = target_name, values_from = mean_ct) %>%
  group_by(sample_name, Task) %>%
  mutate(dCT = DUSP11 - RNA18S1)

# Calculate mock mean dCT of 18S
mock_mean <- data_dct %>%
  filter(sample_name == "mock") %>%
  pull(dCT) %>%
  mean()

# Calculate ddCT by subtracting the dCT of the mock sample from each condition
data_ddct <- data_dct %>%
  group_by(Task) %>%
  mutate(ddCT = dCT - mock_mean) %>%
  ungroup()

# Calculate fold change using the formula 2^(-ddCT)
data_fold_change <- data_ddct %>%
  mutate(fold_change = 2^(-ddCT))

# Calculate welch anova
pvals <- data_fold_change %>%
  # Filter infected samples to remove NA
  filter(sample_name != "infect") %>%
  welch_anova_test(fold_change ~ sample_name)

# Calculate gaines post-hoc test
posthoc <- data_fold_change %>%
  # Again, filter infected samples to remove NA
  filter(sample_name != "infect") %>%
  mutate(sample_name = factor(sample_name, levels = c("mock", "IFN", "polyI:C"))) %>%
  games_howell_test(fold_change ~ sample_name) %>%
  filter(group1 == "mock") %>%
  mutate(y.position = 1.1)

# Export to CSV
write_csv(data_fold_change, "TN072/72.1/72.1_results.csv")

# Plot
plot_data <- data_fold_change %>%
  select(sample_name, Task, fold_change) %>%
  mutate(sample_name = factor(sample_name, levels = c("mock", "IFN", "polyI:C", "infect")))

ggplot(plot_data, aes(x = sample_name, y = fold_change)) +
  geom_point(stat = "identity", position = position_dodge(width=0.5)) +
  labs(title = "DUSP11 mRNA during IFN, polyI:C, and infection", y = "Fold Change") +
  theme_classic() +
  ylim(0, 1.3) +
  theme(legend.title = element_blank(),
  		axis.title.x = element_blank()) +
		stat_pvalue_manual(posthoc, label = "p.adj.signif", hide.ns = FALSE, x = "group2", y.position = 1.3) +
ggsave("TN072/72.1/72.1_plot.png", width = 8, height = 6, dpi = 300)