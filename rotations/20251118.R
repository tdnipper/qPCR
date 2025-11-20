library(tidyverse)
library(dplyr)

df <- read.csv("rotations/2025-11-18_MR_ISG_DUSP11_filtered.csv")

# Calculate mean CT, pivot wider, and convert to numeric
data <- df %>%
  separate(Sample.Name, into = c("Sample.Name", "induction"), sep = "_") %>%
  group_by(Sample.Name, induction, Target.Name, Task) %>%
  summarise(CT = mean(CT, na.rm = TRUE), .groups = "drop")

# Separate induced vs uninduced, calc ddct and fc for each
# during infection (infected vs mock)

data_induced <- data %>%
  filter(induction == "induced") %>%
  pivot_wider(names_from = Target.Name, values_from = CT) %>%
  mutate(across(.cols = c("DUSP11", "ISG15", "MX1", "WSN_PB2"),
                .fns = ~ .x - RNA18S1,
                .names = "delta_ct_{.col}")) %>%
  pivot_longer(
    cols = starts_with("delta_ct_"),
    names_to = "Target",
    values_to = "delta_ct_induced",
    names_prefix = "delta_ct_"
  ) %>%
  select(Sample.Name, induction, Task, Target, delta_ct_induced) %>%
  pivot_wider(names_from = Sample.Name, values_from = delta_ct_induced, names_prefix = "delta_ct_") %>%
  mutate(ddct = delta_ct_infected - delta_ct_mock) %>%
  mutate(fc = 2^(-ddct))

data_uninduced <- data %>%
  filter(induction != "induced") %>%
  pivot_wider(names_from = Target.Name, values_from = CT) %>%
  mutate(across(.cols = c("DUSP11", "ISG15", "MX1", "WSN_PB2"),
                .fns = ~ .x - RNA18S1,
                .names = "delta_ct_{.col}")) %>%
  pivot_longer(
    cols = starts_with("delta_ct_"),
    names_to = "Target",
    values_to = "delta_ct_uninduced",
    names_prefix = "delta_ct_"
  ) %>%
  select(Sample.Name, induction, Task, Target, delta_ct_uninduced) %>%
  pivot_wider(names_from = Sample.Name, values_from = delta_ct_uninduced, names_prefix = "delta_ct_") %>%
  mutate(ddct = delta_ct_infected - delta_ct_mock) %>%
  mutate(fc = 2^(-ddct))

# Combine induced and uninduced data
ddct_infection <- bind_rows(data_induced, data_uninduced) %>%
  select(induction, Task, Target, ddct, fc)

# Separate into mock and infected, calc ddct and fc for each
# during induction (induced vs uninduced)

data_infected <- data %>%
  filter(Sample.Name == "infected") %>%
  pivot_wider(names_from = Target.Name, values_from = CT) %>%
  mutate(across(.cols = c("DUSP11", "ISG15", "MX1", "WSN_PB2"),
                .fns = ~ .x - RNA18S1,
                .names = "delta_ct_{.col}")) %>%
  pivot_longer(
    cols = starts_with("delta_ct_"),
    names_to = "Target",
    values_to = "delta_ct_infected",
    names_prefix = "delta_ct_"
  ) %>%
  select(Sample.Name, induction, Task, Target, delta_ct_infected) %>%
  pivot_wider(names_from = induction, values_from = delta_ct_infected, names_prefix = "delta_ct_") %>%
  mutate(ddct = delta_ct_induced - delta_ct_uninduced) %>%
  mutate(fc = 2^(-ddct))

data_mock <- data %>%
  filter(Sample.Name == "mock") %>%
  pivot_wider(names_from = Target.Name, values_from = CT) %>%
  mutate(across(.cols = c("DUSP11", "ISG15", "MX1", "WSN_PB2"),
                .fns = ~ .x - RNA18S1,
                .names = "delta_ct_{.col}")) %>%
  pivot_longer(
    cols = starts_with("delta_ct_"),
    names_to = "Target",
    values_to = "delta_ct_mock",
    names_prefix = "delta_ct_"
  ) %>%
  select(Sample.Name, induction, Task, Target, delta_ct_mock) %>%
  pivot_wider(names_from = induction, values_from = delta_ct_mock, names_prefix = "delta_ct_") %>%
  mutate(ddct = delta_ct_induced - delta_ct_uninduced) %>%
  mutate(fc = 2^(-ddct))

ddct_induction <- bind_rows(data_infected, data_mock) %>%
  select(Sample.Name, Task, Target, ddct, fc)

# Export to CSV
write.csv(ddct_infection, "rotations/2025-11-18_MR_ddct_infection_results.csv")
write.csv(ddct_induction, "rotations/2025-11-18_MR_ddct_induction_results.csv")

# Plot
p_infection <- ggplot(ddct_infection, aes(x = Target, y = fc, color = induction)) +
  geom_jitter(size = 2, width = 0.05, height = 0.05) +
  labs(title = "Foldchange upon Infection",
       x = "Target Gene",
       y = "Fold Change (18S)") +
  scale_y_log10() +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
print(p_infection)
ggsave("rotations/2025-11-18_MR_ddct_infection_plot.png", plot = p_infection, width = 10, height = 6, dpi = 300)

p_induction <- ggplot(ddct_induction, aes(x = Target, y = fc, color = Sample.Name)) +
  geom_jitter(size = 2, width = 0.05, height = 0.05) +
  labs(title = "Foldchange upon Induction",
       x = "Target Gene",
       y = "Fold Change (18S)") +
  scale_y_log10() +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12))
print(p_induction)
ggsave("rotations/2025-11-18_MR_ddct_induction_plot.png", plot = p_induction, width = 10, height = 6, dpi = 300)



