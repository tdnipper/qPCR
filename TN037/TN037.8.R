library(tidyverse)

# Import data
df <- read.csv("TN037/TN037.8_filtered.csv")

# Summarize mean CT values
data <- df %>%
  group_by(Sample.Name, Target.Name) %>%
    summarize(
        Mean.CT = mean(CT, na.rm = TRUE), 
        .groups = 'drop'
        ) %>%
        extract(Sample.Name, into = c("Sample.Name", "enrichment"), regex = "(.*)\\s(\\S+)$") %>%
        pivot_wider(names_from = enrichment, values_from = Mean.CT)
# Calculate dCT, foldchange, and percent recovery
data <- data %>%
    mutate(
        across(
            .cols = enriched,
            .fns = ~ .x - input  + log2(0.2^-1),
            .names = "dCT_{.col}"
        )
    ) %>%
    mutate(across(
        .cols = starts_with("dCT_"),
        .fns = ~ 2^(-.x),
        .names = "foldchange_{.col}"
    )) %>%
    rename_with(
        .cols = starts_with("foldchange_"),
        .fn = ~ str_replace(.x, "foldchange_dCT_", "foldchange_")
    ) %>%
    mutate(
        across(
            .cols = starts_with("foldchange_"),
            .fns = ~ .x * 100,
            .names = "percent_input_{.col}"
        )) %>%
    rename_with(
        .cols = starts_with("percent_input_"),
        .fn = ~ str_replace(.x, "percent_input_foldchange_", "percent_input_")
    )
write.csv(data, "TN037/TN037.8_recovery_calculations.csv", row.names = FALSE)

# plot percent recovery
plot_data <- data %>%
    select(Target.Name, starts_with("percent_input_"), Sample.Name) %>%
    rename_with(
        .cols = starts_with("percent_input_"),
        .fn = ~ str_replace(.x, "percent_input_", "")
    ) %>%
    pivot_longer(
        cols = enriched,
        names_to = "condition",
        values_to = "percent_input"
    )
p <- ggplot(plot_data, aes(x = Target.Name, y = percent_input, color = Sample.Name)) +
    geom_point(stat = "identity", position = position_dodge(width = 0.5)) +
    labs(
        title = "Percent Input Recovery for TN037.8",
        x = "Target Name",
        y = "Percent Input Recovery",
        color = ""
        ) +
    theme_minimal()
print(p)
ggsave("TN037/TN037.8_percent_input_recovery_plot.png", plot = p, width = 8, height = 6)
