library(tidyverse)

data <- read_csv("TN065/65.8/TN065.8_filtered.csv")

# Take geo mean of technical reps
data <- data |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  summarize(mean_ct = exp(mean(log(CT), na.rm = TRUE)), .groups = "drop") |>
  pivot_wider(names_from = `Target Name`, values_from = mean_ct)

# get geometric mean of PUM1 and TBP for combined housekeeping control
data <- data |>
  rowwise() |>
  mutate(geo_mean = exp(mean(log(c_across(c(PUM1, TBP)))))) |>
  ungroup()

# create dct columns
cols_to_subtract <- c("DUSP11", "WSN_PB2")
data <- data |>
  mutate(across(all_of(cols_to_subtract), 
                list(dct = ~ . - geo_mean),
                .names = "dct {.col}"))
# data <- data |>
#   mutate(across(all_of(cols_to_subtract), 
#                 list(dct = ~ . - PUM1),
#                 .names = "dct {.col}"))

# pivot longer to get dct columns in one column
data_dct <- data |>
  pivot_longer(cols = starts_with("dct"), names_to = "Target", values_to = "dct_value") |>
  mutate(Target = str_remove(Target, "dct ")) |>
  select(`Sample Name`, Task, Target, dct_value)

# get mean dct of DUSP11 and WSN_PB2 for T0 as control
control_dct <- data_dct |>
  filter(`Sample Name` == "T0") |>
  group_by(Target) |>
  summarize(control_dct = mean(dct_value, na.rm = TRUE), .groups = "drop")

# join control dct back to data_dct and calculate ddct
data_ddct <- data_dct |>
  left_join(control_dct, by = "Target") |>
  mutate(ddct = dct_value - control_dct)

# get fold change by calculating 2^-ddct
data_foldchange <- data_ddct |>
  mutate(fold_change = 2^(-ddct))

# save results to csv
write_csv(
  data_foldchange |>
    mutate(across(c(dct_value, control_dct, ddct, fold_change),
    ~ signif(.x, 2)))
  ,"TN065/65.8/foldchange_results.csv"
)

# data to plot foldchange
plot_data <- data_foldchange |>
  select(`Sample Name`, Task, Target, fold_change)

# order x axis by sample name
plot_data$`Sample Name` <- factor(plot_data$`Sample Name`, levels = c("T0", "T8", "T24", "T48"))
plot_data$Task <- as.factor(plot_data$Task)
# take means of foldchange for geom_col later
means <- plot_data |>
  group_by(`Sample Name`, Target) |>
  summarize(mean = mean(fold_change, na.rm = TRUE),
            sd = sd(fold_change, na.rm = TRUE),
            .groups = "drop")

# plot fold change with error bars
p_dusp11 <- ggplot(means %>% filter(Target == "DUSP11"), aes(x = `Sample Name`, y = mean)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), position = position_dodge(width = 0.9), width = 0.2) +
  labs(title = "Fold Change of DUSP11", x = "Sample Name", y = "Fold Change", caption = "Control: TBP+PUM1, 65.8") +
  ylim(0, 1.25) +
  theme_minimal()
ggsave("TN065/65.8/fold_change_DUSP11_TBP_PUM1.png", p_dusp11, width = 6, height = 4, dpi = 300)


amp_data <- read_csv("TN065/65.8/TN065.8_amp.csv") |>
  group_by(`Sample Name`, `Target Name`, Cycle, Task) |>
  summarize(meanRn = mean(`Delta Rn`), sdRN = sd(`Delta Rn`), .groups = "drop")
  
plot_amp <- ggplot(amp_data, aes(x = Cycle, y = meanRn, color = `Sample Name`, group = interaction(`Sample Name`, Task))) +
  geom_line() +
  geom_ribbon(aes(ymin = meanRn - sdRN, ymax = meanRn + sdRN), alpha = 0.2) +
  geom_vline(xintercept=30, linetype="dashed", color = "red") +
  geom_vline(xintercept=4, linetype="dashed", color = "red") +
  geom_rect(aes(xmin=30.1, xmax=Inf, ymin=-Inf, ymax=Inf), fill="grey", color=NA, alpha=0.01) +
  labs(title = "Amplification Curves", x = "Cycle", y = "Mean Delta Rn", caption = "65.8") +
  theme_classic() +
  facet_wrap(~ `Target Name`)
  
ggsave("TN065/65.8/amplification_curves.png", plot_amp, width = 6, height = 4, dpi = 300)