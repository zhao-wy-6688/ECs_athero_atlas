library(ggplot2)
library(dplyr)
library(tidyr)
library(splines)
# 1. Data transformation
score_long <- score %>%
  pivot_longer(cols = c("St", "Ss", "Rts"), names_to = "variable", values_to = "expression")
# 2. Fit smoothing models for "St" and "Ss"
fit_tip <- lm(expression ~ splines::ns(Monocle3_Pseudotime, 3), data = filter(score_long, variable == "St"))
fit_stalk <- lm(expression ~ splines::ns(Monocle3_Pseudotime, 3), data = filter(score_long, variable == "Ss"))
# 3. Generate prediction data
df_pred <- data.frame(Monocle3_Pseudotime = seq(0, 17, length.out = 100))
df_pred$St_smooth <- predict(fit_tip, newdata = df_pred)
df_pred$Ss_smooth <- predict(fit_stalk, newdata = df_pred)
# 4. Identify the crossing point (where the two curves are closest)
df_pred <- df_pred %>%
  mutate(diff = abs(St_smooth - Ss_smooth))  # Calculate absolute difference
cross_point <- df_pred %>%
  filter(diff == min(diff)) %>%
  select(Monocle3_Pseudotime)
# Print the x-coordinate of the crossing point and visualization
ggplot(score_long, aes(x = Monocle3_Pseudotime, y = expression, color = variable)) +
  geom_smooth(method = "lm", formula = y ~ splines::ns(x, 3)) +  
  xlim(0, 17)+
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_manual(values = c("black", "red","darkorange"))+
  geom_vline(xintercept = 5.68, linetype = "dashed", color = "#f18d38", size = 1) + 
  geom_vline(xintercept = 7.35, linetype = "dashed", color = "#c34a44", size = 1) +  
  geom_vline(xintercept = 12.36, linetype = "dashed", color = "#4a6681", size = 1)  ）# cross_point