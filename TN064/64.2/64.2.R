library(tidyverse)

df <- read.csv("TN064/64.2/TN064.2_filtered.csv")

data <- df %>%
  group_by(Sample.Name, Target.Name) %>%
  summarize(mean_Ct = mean(CT, na.rm = TRUE), .groups = "drop") %>%
  separate(Sample.Name, into = c("oligo", "concentration", "condition"), sep = "_")

final <- data %>%
  pivot_wider(names_from = condition, values_from = mean_Ct) %>%
  mutate(delta_ct = enrich - input + log2(0.1^-1)) %>%
  mutate(fold_change = 2^-(delta_ct))

 p <- ggplot(final, aes(x = Target.Name, y = fold_change, color = concentration)) +
   geom_point(aes(shape=oligo)) +
   theme_classic() +
   labs(title = "PB2 Bio-ASO Depletion",
        x = "Transcript",
        y = "Fold Change",
        caption = "64.2"
   )
 print(p)
 ggsave("TN064/64.2/TN064.2_foldchange_plot.png", plot = p, width = 8, height = 6, dpi = 300) 
 
 