library(tidyverse)

# read in data
data <- read_csv("TN069/69.15/TN069.15_filtered.csv")

data <- data |>
  group_by(`Sample Name`, `Target Name`, Task) |>
  summarize(mean_ct = mean(CT), .groups = "drop")
  
# calculate dCT of flowthrough and enrich from input for each sample and target
data_dct <- data |>
  pivot_wider(names_from = Task, values_from = mean_ct) |>
  mutate(dCT_flowthrough = `flowthrough` - `input` + log2(0.1^-1), # input is 10% of total, so adjust by log2(0.1^-1)
         dCT_enrich = `enrich` - `input` + log2(0.1^-1)) # input is 10% of total, so adjust by log2(0.1^-1)

data_recovery <- data_dct |>
  mutate(percent_input_flowthrough = 2^(-dCT_flowthrough) * 100, # convert dCT to percent input
         percent_input_enrich = 2^(-dCT_enrich) * 100) |> # convert dCT to percent input
  select(`Sample Name`, `Target Name`, percent_input_flowthrough, percent_input_enrich)

# save percent input data
write_csv(data_recovery, "TN069/69.15/TN069.15_percent_input.csv")

# plot percent input recovery
plot_data <- data_recovery |>
  pivot_longer(cols = starts_with("percent_input"), names_to = "condition", values_to = "percent_input") |>
  mutate(condition = str_replace(condition, "percent_input_", ""))

p <- ggplot(plot_data, aes(x = `Sample Name`, y = percent_input, color = condition, shape = `Target Name`)) +
  geom_point() +
  labs(title = "Recovery Varying VCE and GTP",
       x = "",
       y = "Percent Input",
       caption = "69.15, 5ug RNA from infected cell with 5pg FLUC spike in") +
  ylim(0,100) +
  theme_classic()

ggsave("TN069/69.15/TN069.15_percent_input_plot.png", plot = p, width = 8, height = 6, dpi = 300)
