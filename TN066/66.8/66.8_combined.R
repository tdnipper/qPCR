library(tidyverse)

df1 <- read.csv("TN066/66.8/TN066.8_filtered.csv")
df2 <- read.csv("TN066/66.8/TN066.8_rep2_filtered.csv")

df <- bind_rows(df1, df2)

data <- df %>%
	group_by(Sample.Name, Target.Name, Task) %>%
	summarise(CT = mean(CT, na.rm = TRUE), .groups = 'drop') %>%
	pivot_wider(names_from = Sample.Name, values_from = CT)

data_enrich <- data %>%
	mutate(mock_FLAG_norm = mock_FLAG - mock_IgG,
		infect_FLAG_norm = infect_FLAG - infect_IgG,
		mock_m2_norm = mock_M2 - mock_IgG,
		infect_m2_norm = infect_M2 - mock_IgG) %>%
	mutate(across(.cols = c(mock_FLAG_norm, infect_FLAG_norm, mock_m2_norm, infect_m2_norm),
		~ 2^-.x))

data_enrich_long <- data_enrich %>%
	select(Target.Name, Task, mock_FLAG_norm, infect_FLAG_norm, mock_m2_norm, infect_m2_norm) %>%
	pivot_longer(cols = c(infect_m2_norm, mock_m2_norm, infect_FLAG_norm, mock_FLAG_norm),
		names_to = "Condition",
		values_to = "normalized_foldchange",
		names_prefix = "foldchange_")

p <- ggplot(data_enrich_long, aes(x = Target.Name, y = normalized_foldchange, color = Condition)) +
	geom_point(size = 2, position = position_dodge(width = 0.5)) +
	theme_classic() +
	ylim(0, 5)
ggsave("TN066/66.8/66.8_combined_foldchange.png", plot = p, width = 10, height = 6, dpi=300)