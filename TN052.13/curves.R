library(readxl)
library(dplyr)
library(writexl)

df <- read_excel("TN052.13/TN052_13_redo.xlsx", sheet = "Amplification Data", skip = 46)
names <- read_excel("TN052.13/TN052_13_redo.xlsx", sheet = "Sample Setup", skip = 46) %>%
  select(`Well`, `Sample Name`, `Custom Task`) %>%
  filter(!is.na(`Sample Name`))

df <- df %>%
  select(`Well`, `Cycle`, `Target Name`, Rn) %>%
  filter(!is.na(`Target Name`), !is.na(Rn))
df <- left_join(df, names, by = "Well")
df_mean <- df %>%
  group_by(`Sample Name`, `Target Name`, Cycle, `Custom Task`) %>%
  summarize(
        Rn_mean = mean(Rn, na.rm = TRUE),
        Rn_sd = sd(Rn, na.rm = TRUE),
        .groups = "drop"
  )

df_melt <- read_excel("TN052.13/TN052_13_redo.xlsx", sheet = "Melt Curve Raw Data", skip = 46) %>%
  select(`Well`, `Reading`, `Temperature`, `Derivative`)
melt_names <- read_excel("TN052.13/TN052_13_redo.xlsx", sheet = "Sample Setup", skip = 46) %>%
  select(`Well`, `Sample Name`, `Target Name`, `Custom Task`)
df_melt <- left_join(df_melt, melt_names, by = "Well")


# Plot
library(ggplot2)
p <- ggplot(df_mean, aes(x = Cycle, y = Rn_mean, color = `Custom Task`)) +
  geom_line() +
  geom_ribbon(aes(ymin = Rn_mean - Rn_sd, ymax = Rn_mean + Rn_sd, fill = `Custom Task`), alpha = 0.2, color = NA, show.legend = FALSE) +
  facet_grid(`Target Name` ~ `Sample Name`) +
  labs(title = "Amplification Curves",
       x = "Cycle",
       y = "Normalized Reporter (Rn)",
    ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
# Save the plot
ggsave("TN052.13/amplification_curves.png", plot = p, width = 12, height = 8, dpi = 300)

# Melt curve plot
p_melt <- ggplot(df_melt, aes(x = Temperature, y = Derivative, color = `Custom Task`)) +
  geom_line() +
  facet_wrap(`Sample Name` ~ `Target Name`) +
  labs(title = "Melt Curves",
       x = "Temperature (°C)",
       y = "Derivative",
       color = "Replicate") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
# Save the melt curve plot
ggsave("TN052.13/melt_curves.png", plot = p_melt, width = 12, height = 8, dpi = 300)