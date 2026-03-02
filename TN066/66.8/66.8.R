library(tidyverse)

df <- read.csv("TN066/66.8/TN066.8_filtered.csv")

data <- df %>%
	group_by(Sample.Name, Target.Name, Task) %>%
	summarise(CT = mean(CT, na.rm = TRUE), .groups = 'drop') %>%
	pivot_wider(names_from = Sample.Name, values_from = CT)	

data_enrich <- data %>%
	mutate(mock_FLAG = mock_FLAG - mock_IgG,
		infect_FLAG = infect_FLAG - infect_IgG,
		mock_m2 = mock_M2 - mock_IgG,
		infect_m2 = infect_M2 - infect_IgG) %>%
	mutate(across(.cols = c(mock_FLAG, infect_FLAG, mock_m2, infect_m2),
		~ 2^-.x))	

data_enrich_long <- data_enrich %>%
    select(Target.Name, Task, mock_FLAG, infect_FLAG, mock_m2, infect_m2) %>%
    pivot_longer(cols = c(infect_m2, mock_m2),
        names_to = "Condition",
        values_to = "normalized_foldchange",
        names_prefix = "foldchange_")

p <- ggplot(data_enrich_long, aes(x = Target.Name, y = normalized_foldchange, color = Condition)) +
	geom_point(size = 2, position = position_dodge(width = 0.5)) +
	theme_classic()
ggsave("TN066/66.8/66.8_foldchange.png", plot = p, width = 10, height = 6, dpi=300)


df_amp <- read.csv("TN066/66.8/TN066.8_rep1_amp.csv") %>%
	select(-Well) %>%
	drop_na() %>%
	group_by(Sample.Name, Target.Name, Cycle, Task) %>%
	summarize(
		Rn_mean = mean(Rn, na.rm = TRUE),
		Rn_sd = sd(Rn, na.rm = TRUE),
		.groups = "drop"
	)
p_amp <- ggplot(df_amp, aes(x = Cycle, y = Rn_mean, color = Sample.Name)) +
	geom_line() +
	geom_ribbon(aes(ymin = Rn_mean - Rn_sd, ymax = Rn_mean + Rn_sd, fill = Sample.Name), alpha = 0.2, color = NA, show.legend = FALSE) +
	theme_classic() + 
	facet_grid(Target.Name ~ Sample.Name)
ggsave("TN066/66.8/66.8_amp.png", plot = p_amp, width = 10, height = 6, dpi=300)
