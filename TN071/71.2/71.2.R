library(tidyverse)

# load in amp data to visualize curves
data <- read_csv("TN071/71.2/TN071.2_amp.csv")

 # take mean and SD of technical replicates
data <- data |>
	group_by(`Sample Name`, `Target Name`, Task, Cycle) |>
	summarize(meanRn = mean(`Delta Rn`), sdRn = sd(`Delta Rn`), .groups = "drop")

# plot amplification curves
data$`Target Name` <- factor(data$`Target Name`, levels = c("RNA18S1", "TBP", "HPRT1", "PUM1", "PPIA", "DUSP11", "WSN_PB2"))
p_amp <- ggplot(data, aes(x = Cycle, y = meanRn, color = `Sample Name`, group = interaction(`Sample Name`, Task))) +
	geom_line() +
	geom_ribbon(aes(ymin = meanRn - sdRn, ymax = meanRn + sdRn, fill = `Sample Name`), alpha = 0.2) +
	labs(
		title = "Amplification Curves",
		subtitle = "Average technical replicates",
		x = "",
		y = "Delta Rn",) +
	theme_classic() +
	facet_wrap(~ `Target Name`)

ggsave("TN071/71.2/TN071.2_amplification_plot.png", plot = p_amp, width = 8, height = 6, dpi = 300)