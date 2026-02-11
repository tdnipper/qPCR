library(tidyverse)

df <- read.csv("TN069/69.8/TN069.8_filtered.csv")

data <- df %>%
    group_by(Sample.Name, Target.Name, Task) %>%
    summarize(
        mean.CT = mean(CT, na.rm = TRUE),
        .groups = 'drop'
    ) %>%
    pivot_wider(names_from = Sample.Name, values_from = mean.CT)

data_enrich <- data %>%
    mutate(
        across(
            .cols = c(enrich, flowthrough),
            .fns = ~ .x - input + log2(0.1^-1),
            .names = "dCT_{.col}"
        )
    ) %>%
    mutate(
        across(
            .cols = starts_with("dCT_"),
            .fns = ~ 2^(-.x),
            .names = "foldchange_{.col}"
        )
    ) %>%
    rename_with(
        .cols = starts_with("foldchange_"),
        .fn = ~ str_replace(.x, "foldchange_dCT", "foldchange")
    ) %>%
    mutate(
        across(
            .cols = starts_with("foldchange_"),
            .fns = ~ .x * 100,
            .names = "percent_input_{.col}"
        )
    ) %>%
    rename_with(
        .cols = starts_with("percent_input"),
        .fn = ~ str_replace(.x, "percent_input_foldchange_", "percent_input_")
    )
write.csv(data_enrich, "TN069/69.8/TN069.8_output.csv", row.names = FALSE)

plot_data <- data_enrich %>%
    select(Target.Name, starts_with("percent_input"), Task) %>%
    rename_with(
        .cols = starts_with("percent_input"),
        .fn = ~ str_replace(.x, "percent_input_", "")
    ) %>%
    pivot_longer(
        cols = c("enrich", "flowthrough"),
        names_to = "condition",
        values_to = "percent_input"
    )

p <- ggplot(plot_data, aes(x = Target.Name, y = percent_input, color = condition)) +
    geom_point(stat = "identity", position = position_dodge(width = 0.5)) +
    labs(
        title = "Percent Input Recovery TN069.8",
        x = "Target Name",
        y = "Percent Input",
        color = ""
        ) +
    theme_minimal()
ggsave("TN069/69.8/TN069.8_percent_input_recovery_plot.png", plot = p, width = 8, height = 6)
