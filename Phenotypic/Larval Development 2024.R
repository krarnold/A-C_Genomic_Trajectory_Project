# Set data_dir to the path of your local data directory
data_dir <- "."

library(ggplot2)
library(data.table)
library(gridExtra)
library(nlme)

se <- function(data) {
  sd(data) / sqrt(length(data))
}

meanDevelopment <- function(df) {
  temp <- c()
  for (i in unique(df$Time)) {
    temp <- c(temp, rep(i, df$totalIndDev[df$Time == i]))
  }
  temp
}

meanDevelopment2 <- function(df) {
  temp <- c()
  for (i in df$scaledTime) {
    temp <- c(temp, rep(i, df$totalIndDev[df$scaledTime == i]))
  }
  temp
}

pupa2024 <- read.csv(file.path(data_dir, "Pupariation_2024.csv"))

pupa2024$totalIndDev <- rowSums(pupa2024[, 7:9])
pupa2024$popCumSum <- rowSums(pupa2024[,4:6])
pupa2024$totalIndDev <- rowSums(pupa2024[, 7:9])
pupa2024$percDev <- 1
pupa2024$percDevA <- 1
pupa2024$percDevB <- 1
pupa2024$percDevC <- 1

for (i in 1:max(pupa2024$Pop)) {
  pupa2024$percDev[pupa2024$Pop == i] <- pupa2024$popCumSum[pupa2024$Pop == i] / max(pupa2024$popCumSum[pupa2024$Pop == i])
}
for (i in 1:max(pupa2024$Pop)) {
  pupa2024$percDevA[pupa2024$Pop == i] <- pupa2024$CumSumA[pupa2024$Pop == i] / max(pupa2024$CumSumA[pupa2024$Pop == i])
}
for (i in 1:max(pupa2024$Pop)) {
  pupa2024$percDevB[pupa2024$Pop == i] <- pupa2024$CumSumB[pupa2024$Pop == i] / max(pupa2024$CumSumB[pupa2024$Pop == i])
}
for (i in 1:max(pupa2024$Pop)) {
  pupa2024$percDevC[pupa2024$Pop == i] <- pupa2024$CumSumC[pupa2024$Pop == i] / max(pupa2024$CumSumC[pupa2024$Pop == i])
}

pupa2024$SimpleSel <- 1
pupa2024$SimpleSel[pupa2024$Treatment == "ACO"] <- "A"
pupa2024$SimpleSel[pupa2024$Treatment == "AO"] <- "A"
pupa2024$SimpleSel[pupa2024$Treatment == "ANCO"] <- "New A"
pupa2024$SimpleSel[pupa2024$Treatment == "NACO"] <- "New A"
pupa2024$SimpleSel[pupa2024$Treatment == "CACO"] <- "New C"
pupa2024$SimpleSel[pupa2024$Treatment == "CAO"] <- "New C"
pupa2024$SimpleSel[pupa2024$Treatment == "CO"] <- "C"
pupa2024$SimpleSel[pupa2024$Treatment == "NCO"] <- "C"

pupa2024$SimpleSel2 <- 1
pupa2024$SimpleSel2[pupa2024$SimpleSel == "A"] <- "Accelerated"
pupa2024$SimpleSel2[pupa2024$SimpleSel == "New A"] <- "Accelerated"
pupa2024$SimpleSel2[pupa2024$SimpleSel == "C"] <- "Control"
pupa2024$SimpleSel2[pupa2024$SimpleSel == "New C"] <- "Control"

pupa2024$OldVsNew <- 1
pupa2024$OldVsNew[pupa2024$SimpleSel == "A"] <- "Ancestral"
pupa2024$OldVsNew[pupa2024$SimpleSel == "New A"] <- "Recent"
pupa2024$OldVsNew[pupa2024$SimpleSel == "C"] <- "Ancestral"
pupa2024$OldVsNew[pupa2024$SimpleSel == "New C"] <- "Recent"

meanRep0 <- 0
for (i in seq(5, 40, 5)) {
  #print(max(pupa2024$Time[pupa2024$Pop == i & pupa2024$percDev < 0.75]))
  #meanRep0 <- meanRep0 + max(pupa2024$Time[pupa2024$Pop == i & pupa2024$percDev < thresh])
  meanRep0 <- meanRep0 + mean(meanDevelopment(pupa2024[pupa2024$Pop == i,]))
}
meanRep0 <- meanRep0 / 8

meanRep1 <- 0
for (i in seq(1, 40, 5)) {
  #meanRep1 <- meanRep1 + max(pupa2024$Time[pupa2024$Pop == i & pupa2024$percDev < thresh])
  meanRep1 <- meanRep1 + mean(meanDevelopment(pupa2024[pupa2024$Pop == i,]))
}
meanRep1 <- meanRep1 / 8

meanRep2 <- 0
for (i in seq(2, 40, 5)) {
  #meanRep2 <- meanRep2 + max(pupa2024$Time[pupa2024$Pop == i & pupa2024$percDev < thresh])
  meanRep2 <- meanRep2 + mean(meanDevelopment(pupa2024[pupa2024$Pop == i,]))
}
meanRep2 <- meanRep2 / 8

meanRep3 <- 0
for (i in seq(3, 40, 5)) {
  #meanRep3 <- meanRep3 + max(pupa2024$Time[pupa2024$Pop == i & pupa2024$percDev < thresh])
  meanRep3 <- meanRep3 + mean(meanDevelopment(pupa2024[pupa2024$Pop == i,]))
}
meanRep3 <- meanRep3 / 8

meanRep4 <- 0
for (i in seq(4, 40, 5)) {
  #meanRep4 <- meanRep4 + max(pupa2024$Time[pupa2024$Pop == i & pupa2024$percDev < thresh])
  meanRep4 <- meanRep4 + mean(meanDevelopment(pupa2024[pupa2024$Pop == i,]))
}
meanRep4 <- meanRep4 / 8

pupa2024$Replicate <- pupa2024$Pop %% 5

pupa2024$scaledTime <- pupa2024$Time
pupa2024$scaledTime[pupa2024$Replicate == 1] <- pupa2024$Time[pupa2024$Replicate == 1] - meanRep1
pupa2024$scaledTime[pupa2024$Replicate == 0] <- pupa2024$Time[pupa2024$Replicate == 0] - meanRep0
pupa2024$scaledTime[pupa2024$Replicate == 2] <- pupa2024$Time[pupa2024$Replicate == 2] - meanRep2
pupa2024$scaledTime[pupa2024$Replicate == 3] <- pupa2024$Time[pupa2024$Replicate == 3] - meanRep3
pupa2024$scaledTime[pupa2024$Replicate == 4] <- pupa2024$Time[pupa2024$Replicate == 4] - meanRep4

meanAC <- 0
for (i in c(11:25, 36:40)) {
  meanAC <- meanAC + mean(meanDevelopment2(pupa2024[pupa2024$Pop == i,]))
}
meanAC <- meanAC / 20

pupa2024$scaledTime2 <- pupa2024$scaledTime - meanAC

pupa2024$Pop2 <- pupa2024$Pop + 60

textSize <- 16
my_theme = theme(
  axis.title.x = element_text(size = textSize, face = "bold"),
  axis.title.y = element_text(size = textSize, face = "bold"),
  axis.text.x = element_text(size = textSize - 4),
  axis.text.y = element_text(size = textSize - 4),
  legend.title = element_text(size = textSize - 4),
  legend.text = element_text(size = textSize - 4),
  legend.background = element_rect(fill="white",
                                   linewidth=0.5, linetype="solid", 
                                   colour ="black"),
  plot.title = element_text(size=textSize, hjust = 0.5, face = "bold"),
  plot.margin = margin(t = 20, r = 5, b = 5, l = 15, unit = "pt"))

cols <- c("A" = "#BB271A", "C" = "#0000F3", "New C" = "forestgreen", "New A" = "#11A0AB")
c <- pupa2024$OldVsNew == "Ancestral"
c <- pupa2024$OldVsNew == "Recent"
c <- pupa2024$Treatment != 7
pup2024_graph <- ggplot(data = pupa2024[c,], mapping = aes(x = scaledTime + max(scaledTime), y = percDev,
                                                           color = as.factor(SimpleSel),
                                                           group = as.factor(Pop)
)) +
  geom_line(linewidth = 0.7, stat = "summary", fun = "mean", alpha = 0.8) +
  ylab("Proportion Pupariated") +
  scale_color_manual(
    values = cols,
    breaks=c('A', 
             'C',
             'New A',
             'New C'),
    labels=c('Founder A', 
             'Founder C', 
             'C>A',
             'A>C')) +
  xlab("Time (Hours)") +
  xlim(0, 100) +
  theme_light() +
  theme(legend.position = c(.85,.5)) +
  labs(color = "Selection",
       linetype = "Evolutionary\nHistory") + 
  my_theme

pup2024_graph

ggsave(pup2024_graph,
       filename = "Pup_Development_2024.tiff",
       device = "tiff",
       width = 7,
       height = 4.5,
       dpi = 300,
       units = "in")