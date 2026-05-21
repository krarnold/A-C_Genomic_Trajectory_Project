library(ggplot2)
library(ggpubr)
library(nlme)

statsTest <- function(c, numSamples, data, newVsOld = FALSE, 
                      year = FALSE, toLog = TRUE, numDays = 3) {
  maxAge <- max(data$Age[c])
  
  print(maxAge)
  
  #intervals <- 1:ceiling(maxAge / 3)
  intervals <- list()
  
  for (i in seq(2, maxAge, numDays)) {
    tempC <- c & (data$Age < i + numDays & data$Age >= i) &
      data$AgeSpecificMort != 0
    
    smallerDf <- data.frame(data$Age[tempC])
    names(smallerDf)[1] <- "Age"
    if (toLog) { smallerDf$AgeSpecificMort <- log(data$AgeSpecificMort[tempC]) }
    else { smallerDf$AgeSpecificMort <- data$AgeSpecificMort[tempC] }
    if (newVsOld) {smallerDf$Selection <- data$SimpleSel[tempC]}
    else if (year) {smallerDf$Selection <- data$SelYear[tempC]}
    else {smallerDf$Selection <- data$Selection[tempC]}
    smallerDf$Population <- data$Population[tempC]
    
    numPopAlive <- length(data$Age[c & (data$Age < i + numDays & data$Age >= i)])
    if (numPopAlive == numSamples * numDays) {
      test <- lme(AgeSpecificMort ~ 
                    Selection, data = smallerDf, 
                  random = ~1|(Population/Age))
      
      #intervals[ceiling((i - 1) / 3)] <- summary(test)$tTable[, metric][2]
      intervals[[ceiling((i - 1) / numDays)]] <- test
    }
  }
  
  intervals
}

extractStats <- function(listOfLME, numDays, numMax = FALSE) {
  if (numMax == FALSE) {
    numMax <- length(listOfLME)
  }
  result <- data.frame(1:numMax)
  colnames(result) <- 'value'
  result$Std.Error
  result$lower <- 1:numMax
  result$higher <- 1:numMax
  result$p.value <- 1:numMax
  result$day <- seq(1, numMax * numDays, numDays) + 13
  
  for (i in 1:numMax) {
    tempSummary <- summary(listOfLME[[i]])$tTable
    value <- tempSummary[2,1]
    std <- tempSummary[2,2]
    mult <- 2.306004
    result$value[i] <- value
    result$Std.Error[i] <- std
    result$lower[i] <- value - std * mult
    result$higher[i] <- value + std * mult
    result$p.value[i] <- tempSummary[2, 5]
  }
  
  threshold <- .05 / numMax
  result$sig <- result$p.value < threshold
  
  result
}

c <- (mortality2023$SimpleSel == "newA" | mortality2023$SimpleSel == "A") & mortality2023$Sex == "F"
test4 <- statsTest(c, 20, mortality2023, "p-value", newVsOld = TRUE, toLog = FALSE, numDays = 1)
result2023NewAA <- extractStats(test4, 1)

c <- (mortality2023$SimpleSel == "newC" | mortality2023$SimpleSel == "C") & mortality2023$Sex == "F"
test4 <- statsTest(c, 20, mortality2023, "p-value", newVsOld = TRUE, toLog = FALSE, numDays = 1)
result2023NewCC <- extractStats(test4, 1)

c <- (mortality2020$SimpleSel == "newA" | mortality2020$SimpleSel == "A") & mortality2020$Sex == "F"
test4 <- statsTest(c, 20, mortality2020, "p-value", newVsOld = TRUE, toLog = FALSE, numDays = 1)
result2020NewAA <- extractStats(test4, 1)

c <- (mortality2020$SimpleSel == "newC" | mortality2020$SimpleSel == "C") & mortality2020$Sex == "F"
test4 <- statsTest(c, 20, mortality2020, "p-value", newVsOld = TRUE, toLog = FALSE, numDays = 1)
result2020NewCC <- extractStats(test4, 1)

c <- (mortality$SimpleSel == "C" | mortality$SimpleSel == "A") & mortality$Sex == "F"
test4 <- statsTest(c, 20, mortality, "p-value", newVsOld = TRUE, toLog = FALSE, numDays = 1)
result2019AC <- -extractStats(test4, 1)

c <- (mortality2020$SimpleSel == "A" | mortality2020$SimpleSel == "C") & mortality2020$Sex == "F"
test4 <- statsTest(c, 20, mortality2020, "p-value", newVsOld = TRUE, toLog = FALSE, numDays = 1)
result2020AC <- -extractStats(test4, 1)

c <- (mortality2023$SimpleSel == "A" | mortality2023$SimpleSel == "C") & mortality2023$Sex == "F"
test4 <- statsTest(c, 20, mortality2023, "p-value", newVsOld = TRUE, toLog = FALSE, numDays = 1)
result2023AC <- -extractStats(test4, 1)

result2020AC$Sel <- "T0"
result2020NewCC$Sel <- "T1"
result2023NewCC$Sel <- "T2"
resultNewC <- rbind(result2020AC, result2020NewCC, result2023NewCC)
resultNewC$ymin <- resultNewC$value - resultNewC$Std.Error
resultNewC$ymax <- resultNewC$value + resultNewC$Std.Error

wocC <- ggplot(resultNewC, mapping = aes(x = abs(day), y= value, color = Sel)) +
  geom_point(size = 3) +
  geom_smooth(linewidth = 1, se = F) +
  theme_light() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.box.background = element_rect(colour = "black")) +
  labs(#title = "Differences in Age Specific Mortality", 
    x = "Age From Egg", 
    y = "Difference in Mean\nAge-Specific Mortality", 
    color = "Legend",
    shape = "Legend") +  
  scale_fill_discrete(breaks=c('T0', 
                               'T1',
                               'T2')) +
  scale_color_discrete(palette = c('#8F00FF', '#0000F3', '#1A85FF'), 
                       breaks=c('T0', 
                                'T1',
                                'T2'),
                       labels=c('Baseline: Founder A- vs. C-types', 
                                '12 Generations of C-type Selection', 
                                '51 Generations of C-type Selection')) +
  my_theme
wocC

result2020AC_forA <- result2020AC
result2020AC_forA$value <- -result2020AC_forA$value
result2020AC_forA$Sel <- "T0"
result2020NewAA$Sel <- "T1"
result2023NewAA$Sel <- "T2"
resultNewA <- rbind(result2020AC_forA, result2020NewAA, result2023NewAA)
resultNewA$ymin <- resultNewA$value - resultNewA$Std.Error
resultNewA$ymax <- resultNewA$value + resultNewA$Std.Error

wocA <- ggplot(resultNewA, mapping = aes(x = abs(day), y= value, color = Sel)) +
  geom_point(size = 3) +
  geom_smooth(linewidth = 1, se = F) +
  theme_light() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.box.background = element_rect(colour = "black")) + 
  #geom_line(color="blue", size = 1, mapping = aes(seq(1, numMax * 3, 3) + 13, value)) +
  scale_fill_discrete(palette = c('purple', '#DC3220', '#BB271A'),
                      breaks=c('T0', 
                               'T1',
                               'T2')) +
  scale_color_discrete(palette = c('#8F00FF', '#BB271A', '#AA33A6'), 
                       breaks=c('T0', 
                                'T1',
                                'T2'),
                       labels=c('Baseline: Founder C- vs. A-types', 
                                '33 Generations of A-type Selection', 
                                '143 Generations of A-type Selection')) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1),
    legend.box.background = element_rect(colour = "black")) +
  xlim(14, 84) +
  labs(#title = "Differences in Age Specific Mortality", 
    x = "Age From Egg", 
    y = "Difference in Mean\nAge-Specific Mortality", 
    color = "Legend",
    shape = "Legend") +
  my_theme
wocA

x <- ggarrange(wocC, wocA, nrow = 2, ncol = 1, align = "v")
ggsave(x, filename = "2026_05_03_WoC_V1.svg", device = "svg", width = 4000, height = 4250, dpi = 300, units = "px")
x <- ggarrange(wocA, wocC, nrow = 2, ncol = 1, align = "v")
ggsave(x, filename = "2026_05_03_WoC_V1_shuffled.svg", device = "svg", width = 4000, height = 4250, dpi = 300, units = "px")

ggsave(x, filename = "2026_05_03_WoC_V1_shuffled.tiff", device = "tiff", width = 4000, height = 4250, dpi = 300, units = "px")