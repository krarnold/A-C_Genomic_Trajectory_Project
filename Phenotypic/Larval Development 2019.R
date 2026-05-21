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

x <- read.csv("Pupariation_2019_CtoA.csv")
x2 <- read.csv("Pupariation_2019_AtoC.csv")
x$totalIndDev <- rowSums(x[, 7:9])
x$popCumSum <- rowSums(x[,4:6])
x$totalIndDev <- rowSums(x[, 7:9])
x$percDev <- 1
x2$percDevA <- 1
x2$percDevB <- 1
x2$percDevC <- 1

for (i in 1:max(x$Pop)) {
  x$percDev[x$Pop == i] <- x$popCumSum[x$Pop == i] / max(x$popCumSum[x$Pop == i])
}
for (i in 1:max(x$Pop)) {
  x$percDevA[x$Pop == i] <- x$CumSumA[x$Pop == i] / max(x$CumSumA[x$Pop == i])
}
for (i in 1:max(x$Pop)) {
  x$percDevB[x$Pop == i] <- x$CumSumB[x$Pop == i] / max(x$CumSumB[x$Pop == i])
}
for (i in 1:max(x$Pop)) {
  x$percDevC[x$Pop == i] <- x$CumSumC[x$Pop == i] / max(x$CumSumC[x$Pop == i])
}

x$SimpleSel <- 1
x$SimpleSel[x$Treatment == "ACO"] <- "A"
x$SimpleSel[x$Treatment == "AO"] <- "A"
x$SimpleSel[x$Treatment == "ANCO"] <- "New A"
x$SimpleSel[x$Treatment == "NACO"] <- "New A"
x$SimpleSel[x$Treatment == "CACO"] <- "New C"
x$SimpleSel[x$Treatment == "CAO"] <- "New C"
x$SimpleSel[x$Treatment == "CO"] <- "C"
x$SimpleSel[x$Treatment == "NCO"] <- "C"



x2$totalIndDev <- rowSums(x2[, 7:9])
x2$popCumSum <- rowSums(x2[,4:6])
x2$percDev <- 1
x2$percDevA <- 1
x2$percDevB <- 1
x2$percDevC <- 1

for (i in 1:max(x2$Pop)) {
  x2$percDev[x2$Pop == i] <- x2$popCumSum[x2$Pop == i] / max(x2$popCumSum[x2$Pop == i])
}
for (i in 1:max(x2$Pop)) {
  x2$percDevA[x2$Pop == i] <- x2$CumSumA[x2$Pop == i] / max(x2$CumSumA[x2$Pop == i])
}
for (i in 1:max(x2$Pop)) {
  x2$percDevB[x2$Pop == i] <- x2$CumSumB[x2$Pop == i] / max(x2$CumSumB[x2$Pop == i])
}
for (i in 1:max(x2$Pop)) {
  x2$percDevC[x2$Pop == i] <- x2$CumSumC[x2$Pop == i] / max(x2$CumSumC[x2$Pop == i])
}
x2$SimpleSel <- 1
x2$SimpleSel[x2$Treatment == "ACO"] <- "A"
x2$SimpleSel[x2$Treatment == "AO"] <- "A"
x2$SimpleSel[x2$Treatment == "ANCO"] <- "New A"
x2$SimpleSel[x2$Treatment == "NACO"] <- "New A"
x2$SimpleSel[x2$Treatment == "CACO"] <- "New C"
x2$SimpleSel[x2$Treatment == "CAO"] <- "New C"
x2$SimpleSel[x2$Treatment == "CO"] <- "C"
x2$SimpleSel[x2$Treatment == "NCO"] <- "C"

x2$Replicate <- x2$Pop %% 5

meanRep0 <- 0
for (i in seq(5, 30, 5)) {
  #print(max(x2$Time[x2$Pop == i & x2$percDev < 0.75]))
  #meanRep0 <- meanRep0 + max(x2$Time[x2$Pop == i & x2$percDev < thresh])
  meanRep0 <- meanRep0 + mean(meanDevelopment(x2[x2$Pop == i,]))
}
meanRep0 <- meanRep0 / 8

meanRep1 <- 0
for (i in seq(1, 30, 5)) {
  #meanRep1 <- meanRep1 + max(x2$Time[x2$Pop == i & x2$percDev < thresh])
  meanRep1 <- meanRep1 + mean(meanDevelopment(x2[x2$Pop == i,]))
}
meanRep1 <- meanRep1 / 8

meanRep2 <- 0
for (i in seq(2, 30, 5)) {
  #meanRep2 <- meanRep2 + max(x2$Time[x2$Pop == i & x2$percDev < thresh])
  meanRep2 <- meanRep2 + mean(meanDevelopment(x2[x2$Pop == i,]))
}
meanRep2 <- meanRep2 / 8

meanRep3 <- 0
for (i in seq(3, 30, 5)) {
  #meanRep3 <- meanRep3 + max(x2$Time[x2$Pop == i & x2$percDev < thresh])
  meanRep3 <- meanRep3 + mean(meanDevelopment(x2[x2$Pop == i,]))
}
meanRep3 <- meanRep3 / 8

meanRep4 <- 0
for (i in seq(4, 30, 5)) {
  #meanRep4 <- meanRep4 + max(x2$Time[x2$Pop == i & x2$percDev < thresh])
  meanRep4 <- meanRep4 + mean(meanDevelopment(x2[x2$Pop == i,]))
}
meanRep4 <- meanRep4 / 8

x2$scaledTime <- x2$Time
x2$scaledTime[x2$Replicate == 1] <- x2$Time[x2$Replicate == 1] - meanRep1
x2$scaledTime[x2$Replicate == 0] <- x2$Time[x2$Replicate == 0] - meanRep0
x2$scaledTime[x2$Replicate == 2] <- x2$Time[x2$Replicate == 2] - meanRep2
x2$scaledTime[x2$Replicate == 3] <- x2$Time[x2$Replicate == 3] - meanRep3
x2$scaledTime[x2$Replicate == 4] <- x2$Time[x2$Replicate == 4] - meanRep4

meanAC <- 0
for (i in c(1:10, 21:30)) {
  meanAC <- meanAC + mean(meanDevelopment2(x2[x2$Pop == i,]))
}
meanAC <- meanAC / 20

x2$scaledTime2 <- x2$scaledTime - meanAC

x$Replicate <- x$Pop %% 5

meanRep0 <- 0
for (i in seq(5, 30, 5)) {
  #print(max(x$Time[x$Pop == i & x$percDev < 0.75]))
  #meanRep0 <- meanRep0 + max(x$Time[x$Pop == i & x$percDev < thresh])
  meanRep0 <- meanRep0 + mean(meanDevelopment(x[x$Pop == i,]))
}
meanRep0 <- meanRep0 / 8

meanRep1 <- 0
for (i in seq(1, 30, 5)) {
  #meanRep1 <- meanRep1 + max(x$Time[x$Pop == i & x$percDev < thresh])
  meanRep1 <- meanRep1 + mean(meanDevelopment(x[x$Pop == i,]))
}
meanRep1 <- meanRep1 / 8

meanRep2 <- 0
for (i in seq(2, 30, 5)) {
  #meanRep2 <- meanRep2 + max(x$Time[x$Pop == i & x$percDev < thresh])
  meanRep2 <- meanRep2 + mean(meanDevelopment(x[x$Pop == i,]))
}
meanRep2 <- meanRep2 / 8

meanRep3 <- 0
for (i in seq(3, 30, 5)) {
  #meanRep3 <- meanRep3 + max(x$Time[x$Pop == i & x$percDev < thresh])
  meanRep3 <- meanRep3 + mean(meanDevelopment(x[x$Pop == i,]))
}
meanRep3 <- meanRep3 / 8

meanRep4 <- 0
for (i in seq(4, 30, 5)) {
  #meanRep4 <- meanRep4 + max(x$Time[x$Pop == i & x$percDev < thresh])
  meanRep4 <- meanRep4 + mean(meanDevelopment(x[x$Pop == i,]))
}
meanRep4 <- meanRep4 / 8

x$scaledTime <- x$Time
x$scaledTime[x$Replicate == 1] <- x$Time[x$Replicate == 1] - meanRep1
x$scaledTime[x$Replicate == 0] <- x$Time[x$Replicate == 0] - meanRep0
x$scaledTime[x$Replicate == 2] <- x$Time[x$Replicate == 2] - meanRep2
x$scaledTime[x$Replicate == 3] <- x$Time[x$Replicate == 3] - meanRep3
x$scaledTime[x$Replicate == 4] <- x$Time[x$Replicate == 4] - meanRep4

meanAC <- 0
for (i in c(1:5, 11:15, 21:30)) {
  meanAC <- meanAC + mean(meanDevelopment2(x[x$Pop == i,]))
}
meanAC <- meanAC / 20

x$scaledTime2 <- x$scaledTime - meanAC

textSize <- 16
my_theme = theme(
  axis.title.x = element_text(size = textSize, face = "bold"),
  axis.title.y = element_text(size = textSize, face = "bold"),
  axis.text.x = element_text(size = textSize - 4),
  axis.text.y = element_text(size = textSize - 4),
  legend.title = element_text(size = textSize - 4),
  legend.text = element_text(size = textSize - 6),
  legend.background = element_rect(fill="white",
                                   linewidth=0.5, linetype="solid", 
                                   colour ="black"),
  plot.title = element_text(size=textSize, hjust = 0.5, face = "bold"))

cols <- c("A" = "#BB271A", "C" = "#0000F3", "New C" = "forestgreen", "New A" = "#11A0AB")

x$Sel2 <- "NewA"
x2$Sel2 <- "NewC"
pup2019 <- rbind(x, x2)
pup2019$Pop2 <- pup2019$Pop
pup2019$Pop2[pup2019$Sel2 == "NewC"] <- pup2019$Pop2[pup2019$Sel2 == "NewC"] + 30

# Shared plotting function to avoid repetition
make_pup_graph <- function(data) {
  ggplot(data = data, mapping = aes(x = scaledTime - min(scaledTime[totalIndDev > 0]), 
                                    y = percDev,
                                    color = as.factor(SimpleSel),
                                    group = as.factor(Pop))) +
    geom_line(alpha = 0.7, linewidth = 0.7, stat = "summary", fun = "mean") +
    ylab("Proportion Pupariated") +
    xlab("Time (Hours)") +
    xlim(0, 100) +
    theme_light() +
    theme(legend.position = c(.85, .5),
          legend.box.background = element_rect(colour = "black")) +
    scale_color_manual(
      values = cols,
      breaks = c('A', 'C', 'New A', 'New C'),
      labels = c('Founder A', 'Founder C', 'C>A', 'A>C')) +
    labs(color = "Selection") +
    my_theme
}

# NewA graph (C to A selection)
pup2019_graph_NewA <- make_pup_graph(pup2019[pup2019$Sel2 == "NewA", ])
pup2019_graph_NewA

ggsave(pup2019_graph_NewA,
       filename = "Pup_Development2019_NewA.tiff",
       device = "tiff",
       width = 7,
       height = 4.5,
       dpi = 300,
       units = "in")

# NewC graph (A to C selection)
pup2019_graph_NewC <- make_pup_graph(pup2019[pup2019$Sel2 == "NewC", ])
pup2019_graph_NewC

ggsave(pup2019_graph_NewC,
       filename = "Pup_Development2019_NewC.tiff",
       device = "tiff",
       width = 7,
       height = 4.5,
       dpi = 300,
       units = "in")