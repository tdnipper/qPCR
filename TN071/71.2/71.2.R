library(tidyverse)
library(ggh4x)
library(ggpubr)
library(rstatix)

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
    y = "Delta Rn"
  ) +
  theme_classic() +
  facet_wrap(~ `Target Name`)

ggsave("TN071/71.2/TN071.2_amplification_plot.png", plot = p_amp, width = 8, height = 6, dpi = 300)

# read in CT data
ct_data <- read_csv("TN071/71.2/TN071.2_filtered.csv")

# summarize CT data to get mean of technical replicates
data_dCT <- ct_data |>
  group_by(`Sample Name`, Task, `Target Name`) |>
  summarize(meanCT = mean(CT), .groups = "drop") |>
  pivot_wider(
    names_from  = `Target Name`,
    values_from = meanCT
  ) |>
  mutate(
    across(
      .cols = c(DUSP11, WSN_PB2),
      .fns = ~ .x - RNA18S1,
      .names = "dCT_{col}"
    )
  ) |>
  select(`Sample Name`, Task, starts_with("dCT_"))

# calculate mean control dct
control_avg <- data_dCT |>
  filter(`Sample Name` == "mock") |>
  summarize(
    control_avg_DUSP11 = mean(dCT_DUSP11, na.rm = TRUE),
    control_avg_WSN_PB2 = mean(dCT_WSN_PB2, na.rm = TRUE),
    .groups = "drop"
  )

# join control average to data_dCT and calculate dCT and fold change
data_foldchange <- data_dCT |>
  cross_join(control_avg) |>
  mutate(
    ddct_DUSP11 = dCT_DUSP11 - control_avg_DUSP11,
    ddct_WSN_PB2 = dCT_WSN_PB2 - control_avg_WSN_PB2
  ) |>
  mutate(
    foldchange_DUSP11 = 2^(-ddct_DUSP11),
    foldchange_WSN_PB2 = 2^(-ddct_WSN_PB2)
  ) |>
  select(`Sample Name`, Task, starts_with("foldchange_"), starts_with("ddct_"))

# calculate significance of ddCT
data_ddct <- data_foldchange |>
  select(`Sample Name`, Task, starts_with("ddct_")) |>
  pivot_longer(
    cols = starts_with("ddct_"),
    names_to = "Target",
    values_to = "ddct",
    names_prefix = "ddct_"
  )
  
# assumption checks on DUSP11 ddCT (n = 3/group, so these are only indicative)
# rename to a syntactic column — rstatix formula parsing dislikes the space
dusp <- data_ddct |>
  filter(Target == "DUSP11") |>
  rename(group = `Sample Name`)
print(shapiro.test(residuals(aov(ddct ~ group, data = dusp))))  # normality
print(levene_test(dusp, ddct ~ group)) # equal variance

# QQ diagnostic for DUSP11 ddCT residuals
resid_dusp <- residuals(aov(ddct ~ group, data = dusp))
ggqqplot(resid_dusp) +
  labs(
    title = "DUSP11 ddCT — residual QQ",
    caption = "Shapiro-Wilk p = 0.92, Levene p = 0.99"
  )
ggsave("TN071/71.2/TN071.2_DUSP11_qq.png", width = 4, height = 4, dpi = 300)

# significance of ddCT — Dunnett vs mock (PE+ and PE- each vs mock)
set.seed(42)  # Dunnett's p-value adjustment is Monte-Carlo
dunnett_DUSP11 <- dusp |>
  dunnett_test(ddct ~ group, ref.group = "mock") |>
  mutate(Target = "DUSP11", y.position = c(1.05, 1.15))
print(dunnett_DUSP11)

# write to csv
write_csv(
    mutate(
      data_foldchange,
      across(
        .cols = c(foldchange_DUSP11, foldchange_WSN_PB2, ddct_DUSP11, ddct_WSN_PB2),
        .fns = ~ signif(.x, digits = 2)
      )
    ),
    "TN071/71.2/TN071.2_foldchange.csv"
)

# get data for plotting
data_plot <- data_foldchange |>
  pivot_longer(
    cols = starts_with("foldchange_"),
    names_to = "Target",
    values_to = "foldchange",
    names_prefix = "foldchange_"
  )

plot_means = data_plot |>
  group_by(`Sample Name`, Target) |>
  summarize(
    meanfc = mean(foldchange),
    sd = sd(foldchange),
    .groups = "drop"
  )

ggplot(data = plot_means, aes(x = `Sample Name`, y = meanfc)) +
  geom_col(alpha = 0.7, fill = "#2A3752") +
  geom_errorbar(
    aes(ymin = meanfc - sd, ymax = meanfc + sd),
    width = 0.2
  ) +
  geom_jitter(data = data_plot, aes(x = `Sample Name`, y = foldchange), color = "#2A3752", width = 0.1) +
  stat_pvalue_manual(dunnett_DUSP11, x = "group2", label = "p.adj.signif", y.position = 1.25) +
  labs(
    title = "ALI infection",
    x = "Cell Status",
    y = "Fold Change"
  ) +
  theme_classic() +
  facet_wrap(~ Target, scales = "free_y") +
  # ggh4x allows for facetted scales
  facetted_pos_scales(
    y = list(
      Target == "DUSP11" ~ scale_y_continuous(
        limits = c(0, 1.25),
        breaks = seq(0, 1.25, 0.25)
      ),
      Target == "WSN_PB2" ~ scale_y_continuous()
    )
  )

ggsave("TN071/71.2/TN071.2_foldchange_18S_plot.png", width = 8, height = 6, dpi = 300)