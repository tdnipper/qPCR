library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)

df <- read.csv("TN066//66.2/66.2_filtered.csv")

# Calculate mean CT, pivot wider, and convert to numeric
data <- df %>%
    group_by(Sample.Name, Target.Name, Task) %>%
    summarise(CT = mean(CT, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Sample.Name, values_from = CT)

# Calculate enrichment for IgG and NP samples
data_enrich <- data %>%
    mutate(across(.cols = infect_enrich_NP,
                  .fns  = ~ . - infect_input_NP - log2(0.044^-1),
                  .names = "infected_NP_enrichment")) %>%
    mutate(across(.cols = mock_enrich_NP,
                  .fns  = ~ . - mock_input_NP - log2(0.044^-1),
                  .names = "mock_NP_enrichment")) %>%
    mutate(across(.cols = infect_enrich_IgG,
                  .fns  = ~ . - infect_input_IgG - log2(0.044^-1),
                  .names = "infected_IgG_enrichment")) %>%
    mutate(across(.cols = mock_enrich_IgG,
                  .fns  = ~ . - mock_input_IgG - log2(0.044^-1),
                  .names = "mock_IgG_enrichment")) %>%
    mutate(across(.cols = ends_with("enrichment"),
                  .fns  = ~ 2^-.x,
                  .names = "{.col}_percent_input"))
enrich_out <- data_enrich %>%
    select(Target.Name, Task,
           ends_with("percent_input")) %>%
    rename_with(~ sub("_enrichment_percent_input", "", .x), ends_with("percent_input")) %>%
    pivot_longer(
        cols = -c(Target.Name, Task),
        names_to = "Sample_Type",
        values_to = "percent_input",
        names_prefix = "percent_input_"
    )
write.csv(enrich_out, "TN066/66.2/66.2_percent_input_results.csv", row.names = FALSE)
# Separate mock and infected for plotting
enrich_out_mock <- enrich_out %>%
    filter(grepl("mock", Sample_Type))
enrich_out_infect <- enrich_out %>%
    filter(grepl("infect", Sample_Type))

p1_mock <- ggplot(enrich_out_mock, aes(x = Target.Name, y = percent_input, color = Sample_Type)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    labs(title = "Mock Enrichment",
         x = "Target Gene",
         y = "Percent Input",
         color = "") +
    scale_y_log10(labels = label_number()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
ggsave("TN066/66.2/66.2_percent_input_mock.png", plot = p1_mock, width = 10, height = 6, dpi = 300)

p1_infect <- ggplot(enrich_out_infect, aes(x = Target.Name, y = percent_input, color = Sample_Type)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    labs(title = "Infected Enrichment",
         x = "Target Gene",
         y = "Percent Input",
         color = "") +
    scale_y_log10(labels = label_number()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
ggsave("TN066/66.2/66.2_percent_input_infected.png", plot = p1_infect, width = 10, height = 6, dpi = 300)

# Calculate foldchange for mock and infected samples relative to IgG
data_foldchange <- data_enrich %>%
    mutate(infected = infect_enrich_NP - infect_enrich_IgG,
           mock = mock_enrich_NP - mock_enrich_IgG) %>%
    select(Target.Name, Task, infected, mock) %>%
    mutate(across(c(infected, mock), ~ 2^-.x))

# Pivot foldchange results longer for output
data_foldchange_long <- data_foldchange %>%
    select(Target.Name, infected, mock, Task) %>%
    pivot_longer(cols = c(infected, mock),
                 names_to = "Condition",
                 values_to = "relative_percent_input",
                 names_prefix = "percent_input_")
write.csv(data_foldchange_long, "TN066/66.2/66.2_foldchange_results.csv")

p2 <- ggplot(data_foldchange_long, aes(x = Target.Name, y = relative_percent_input, color = Condition)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    labs(title = "Enrichment",
         x = "Target Gene",
         y = "Foldchange / IgG") +
    scale_y_log10(labels = label_number()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
ggsave("TN066/66.2/66.2_foldchange_plot.png", plot = p2, width = 10, height = 6, dpi = 300)

# Prepare publication figure without PB2
pub_data_mock <- enrich_out_mock %>%
    filter(Target.Name != "WSN_PB2")

pub_data_infect <- enrich_out_infect %>%
    filter(Target.Name != "WSN_PB2")

pub_mock <- ggplot(pub_data_mock, aes(x = Target.Name, y = percent_input, color = Sample_Type)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    labs(title = "Mock Enrichment",
         x = "",
         y = "Percent Input",
         color = "") +
    scale_y_log10(labels = label_number()) +
    coord_cartesian(ylim = c(0.01, 100)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
pub_infect <- ggplot(pub_data_infect, aes(x = Target.Name, y = percent_input, color = Sample_Type)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    labs(title = "Infected Enrichment",
         x = "",
         y = "",
         color = "") +
    scale_y_log10(labels = label_number()) +
    coord_cartesian(ylim = c(0.01, 100)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
pub_figure <- ggarrange(pub_mock, pub_infect, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
print(pub_figure)
ggsave("TN066/66.2/66.2_rip_pub_figure.png", plot = pub_figure, width = 12, height = 6, dpi = 300)
