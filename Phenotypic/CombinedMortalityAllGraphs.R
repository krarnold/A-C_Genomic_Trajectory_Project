mortCombined3Years <- rbind(mortality, mortality2020, mortality2023)

textSize <- 24
my_theme = theme(
  axis.title.x = element_text(size = textSize, face = "bold"),
  axis.title.y = element_text(size = textSize, face = "bold"),
  axis.text.x = element_text(size = textSize - 4),
  axis.text.y = element_text(size = textSize - 4),
  legend.title = element_text(size = textSize - 4),
  legend.text = element_text(size = textSize - 6),
  #legend.position = c(.85,.4),
  #legend.box.background = element_rect(colour = "black")), 
  legend.background = element_rect(fill="white",
                                   linewidth=0.5, linetype="solid", 
                                   colour ="black"),
  plot.title = element_text(size=textSize, hjust = 0.5, face = "bold"))


c <- mortCombined3Years$Sex == "F" & 
  #(mortCombined3Years$SimpleSel == "C" | mortCombined3Years$SimpleSel == "A") &
  (mortCombined3Years$AgeSpecificMort != 0 & mortCombined3Years$Alive > 49) &
  (mortCombined3Years$SelYear == "C2019" | mortCombined3Years$SelYear == "A2019")
ggplot(mapping = aes(x = mortCombined3Years$Age[c],
                     y = (mortCombined3Years$AgeSpecificMort[c]),
                     color = mortCombined3Years$SelYear[c])) + 
  theme_light() +
  geom_point(size = 1.5, alpha = .7) + 
  theme(legend.position.inside = c(.85,.5),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "2019 Mortality Assay", 
       x = "Age From Egg", 
       y = "Age Specific Mortality", 
       color = "Legend",
       shape = "Legend") +
  my_theme

ggplot(mapping = aes(x = mortCombined3Years$Age[c] + 13,
                     y = log(mortCombined3Years$AgeSpecificMort[c]),
                     color = mortCombined3Years$SelYear[c])) + 
  theme_light() +
  geom_point(size = 1.5, alpha = .7) + 
  theme(legend.position = c(.85,.5),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "2019 Mortality Assay", 
       x = "Age From Egg", 
       y = "log(Age Specific Mortality)", 
       color = "Legend",
       shape = "Legend") +
  my_theme


c <- mortCombined3Years$Sex == "F" & 
  #(mortCombined3Years$SimpleSel == "C" | mortCombined3Years$SimpleSel == "A") &
  (mortCombined3Years$AgeSpecificMort != 0 & mortCombined3Years$Alive > 49) & 
  ((mortCombined3Years$SelYear == "A2020" | mortCombined3Years$SelYear == "New A 2020") | 
  (mortCombined3Years$SelYear == "C2020" | mortCombined3Years$SelYear == "New C 2020"))

c <- mortCombined3Years$Sex == "F" & 
  #(mortCombined3Years$SimpleSel == "C" | mortCombined3Years$SimpleSel == "A") &
  (mortCombined3Years$AgeSpecificMort != 0 & mortCombined3Years$Alive > 49) & 
  ((mortCombined3Years$SelYear == "A2020" | mortCombined3Years$SelYear == "New A 2020")) #| 
     #(mortCombined3Years$SelYear == "C2020" | mortCombined3Years$SelYear == "New C 2020"))

ggplot(mapping = aes(x = mortCombined3Years$Age[c] + 13,
                     y = log(mortCombined3Years$AgeSpecificMort[c]),
                     color = mortCombined3Years$SelYear[c])) + 
  theme_light() +
  geom_point(size = 3, alpha = .7) + 
  theme(legend.position = c(.85,.25),
        legend.box.background = element_rect(colour = "black")) +
  labs(#title = "Mortality Assay at ~28 Generations", 
       x = "Age From Egg", 
       y = "log(Age Specific Mortality)", 
       color = "Legend",
       shape = "Legend") +
  my_theme

ggplot(mapping = aes(x = mortCombined3Years$Age[c] + 13,
                     y = log(mortCombined3Years$AgeSpecificMort[c]),
                     color = mortCombined3Years$SelYear[c])) + 
  theme_light() +
  geom_point(size = 3, alpha = .7) + 
  theme(legend.position = c(.85,.25),
        legend.box.background = element_rect(colour = "black")) +
  labs(#title = "2020 Mortality Assay", 
       x = "Age From Egg", 
       y = "log(Age Specific Mortality)", 
       color = "Legend",
       shape = "Legend") +
  my_theme


c <- mortCombined3Years$Sex == "F" & 
  #(mortCombined3Years$SimpleSel == "C" | mortCombined3Years$SimpleSel == "A") &
  (mortCombined3Years$AgeSpecificMort != 0 & mortCombined3Years$Alive > 49) & 
  ((mortCombined3Years$SelYear == "New C 2023" | mortCombined3Years$SelYear == "New A 2023") | 
     (mortCombined3Years$SelYear == "C2023" | mortCombined3Years$SelYear == "A2023"))
ggplot(mapping = aes(x = mortCombined3Years$Age[c] + 13,
                     y = (mortCombined3Years$AgeSpecificMort[c]),
                     color = mortCombined3Years$SelYear[c])) + 
  theme_light() +
  geom_point(size = 1.5, alpha = .7) + 
  theme(legend.position = c(.85,.5),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "2023 Mortality Assay", 
       x = "Age From Egg", 
       y = "Age Specific Mortality", 
       color = "Legend",
       shape = "Legend") +
  my_theme

ggplot(mapping = aes(x = mortCombined3Years$Age[c] + 13,
                     y = log(mortCombined3Years$AgeSpecificMort[c]),
                     color = mortCombined3Years$SelYear[c])) + 
  theme_light() +
  geom_point(size = 1.5, alpha = .7) + 
  theme(legend.position = c(.85,.25),
        legend.box.background = element_rect(colour = "black")) +
  labs(title = "2023 Mortality Assay", 
       x = "Age From Egg", 
       y = "log(Age Specific Mortality)", 
       color = "Legend",
       shape = "Legend") +
  my_theme


c <- (mortCombined3Years$SelYear == "A2019" |
  mortCombined3Years$SelYear == "C2019" |
  mortCombined3Years$SelYear == "A2020" |
  mortCombined3Years$SelYear == "C2020" |
  mortCombined3Years$SelYear == "A2023" |
  mortCombined3Years$SelYear == "C2023") &
  mortCombined3Years$AgeSpecificMort != 0 &
  mortCombined3Years$Alive > 50
g <- ggplot(mapping = aes(x = mortCombined3Years$Age[c] + 13,
                     y = log(mortCombined3Years$AgeSpecificMort[c]),
                     color = mortCombined3Years$SelYear[c])) + 
  theme_light() +
  geom_point(size = 1.5, alpha = .7) + 
  theme(legend.position = c(.85,.25),
        legend.box.background = element_rect(colour = "black")) +
  labs(#title = "2023 Mortality Assay", 
       x = "Age From Egg", 
       y = "log(Age Specific Mortality)", 
       color = "Legend",
       shape = "Legend") +
  my_theme

ggsave(g, filename = "Ancestral_Reproducibility_Mortality_v2.png", device = "png", width = 4000, height = 2250, dpi = 300, units = "px")
png("x.png", width = 4000, height = 2250, res = 300)