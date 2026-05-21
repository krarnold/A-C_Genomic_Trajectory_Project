library(nlme)

statsTest <- function(c, numSamples, data, metric = "depricated", 
                      newVsOld = FALSE, year = FALSE, toLog = TRUE, numDays = 3) {
  maxAge <- max(data$Age[c])
  print(maxAge)
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

numDays <- 3

c <- (mortCombined3Years$SelYear == "C2019" | mortCombined3Years$SelYear == "C2023") & mortCombined3Years$Sex == "F"
test4 <- statsTest(c, 20, mortCombined3Years, year = TRUE)
result4 <- extractStats(test4, 3)
write.csv(result4, "C2020vsC2019.csv")

c <- (mortCombined3Years$SelYear == "C2020" | mortCombined3Years$SelYear == "C2023") & mortCombined3Years$Sex == "F"
test5 <- statsTest(c, 20, mortCombined3Years, year = TRUE)
result5 <- extractStats(test5, 3)
write.csv(result5, "C2020vsC2023.csv")

c <- (mortCombined3Years$SelYear == "C2019" | mortCombined3Years$SelYear == "C2020") & mortCombined3Years$Sex == "F"
test6 <- statsTest(c, 20, mortCombined3Years, year = TRUE)
result6 <- extractStats(test6, 3)
write.csv(result6, "C2019vsC2023.csv")

c <- (mortCombined3Years$SelYear == "A2019" | mortCombined3Years$SelYear == "A2023") & mortCombined3Years$Sex == "F"
test4 <- statsTest(c, 20, mortCombined3Years, year = TRUE)
result4 <- extractStats(test4, 3)
write.csv(result4, "A2020vsA2019.csv")

c <- (mortCombined3Years$SelYear == "A2020" | mortCombined3Years$SelYear == "A2023") & mortCombined3Years$Sex == "F"
test5 <- statsTest(c, 20, mortCombined3Years, year = TRUE)
result5 <- extractStats(test5, 3)
write.csv(result5, "A2020vsA2023.csv")

c <- (mortCombined3Years$SelYear == "A2019" | mortCombined3Years$SelYear == "A2020") & mortCombined3Years$Sex == "F"
test6 <- statsTest(c, 20, mortCombined3Years, year = TRUE)
result6 <- extractStats(test6, 3)
write.csv(result6, "A2019vsA2023.csv")

c <- (mortality2020$SimpleSel == "newC" | mortality2020$SimpleSel == "C") & mortality2020$Sex == "F"
test4 <- statsTest(c, 20, mortality2020, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 1)
result0 <- extractStats(test4, 1)
write.csv(result0, "NewC2020vsC2020.csv")

c <- (mortality2020$SimpleSel == "C" | mortality2020$SimpleSel == "A") & mortality2020$Sex == "F"
test4 <- statsTest(c, 20, mortality2020, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result1 <- extractStats(test4, 3)
write.csv(result1, "C2020vsA2020.csv")

c <- (mortality2020$SimpleSel == "newC" | mortality2020$SimpleSel == "A") & mortality2020$Sex == "F"
test4 <- statsTest(c, 20, mortality2020, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "NewC2020vsA2020.csv")

c <- (mortality2020$SimpleSel == "newA" | mortality2020$SimpleSel == "A") & mortality2020$Sex == "F"
test4 <- statsTest(c, 20, mortality2020, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "NewA2020vsA2020.csv")

c <- (mortality2020$SimpleSel == "newA" | mortality2020$SimpleSel == "C") & mortality2020$Sex == "F"
test4 <- statsTest(c, 20, mortality2020, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "NewA2020vsC2020.csv")

c <- (mortality2023$SimpleSel == "newA" | mortality2023$SimpleSel == "A") & mortality2023$Sex == "F"
test4 <- statsTest(c, 20, mortality2023, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "NewA2023vsA2023.csv")

c <- (mortality2023$SimpleSel == "newC" | mortality2023$SimpleSel == "C") & mortality2023$Sex == "F"
test4 <- statsTest(c, 20, mortality2023, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "NewC2023vsC2023.csv")

c <- (mortality$SimpleSel == "C" | mortality$SimpleSel == "A") & mortality$Sex == "F"
test4 <- statsTest(c, 20, mortality, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "C2019vsA2019.csv")

c <- (mortality2023$SimpleSel == "A" | mortality2023$SimpleSel == "C") & mortality2023$Sex == "F"
test4 <- statsTest(c, 20, mortality2023, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "C2023vsA2023.csv")

c <- (mortality2023$SimpleSel == "newC" | mortality2023$SimpleSel == "A") & mortality2023$Sex == "F"
test4 <- statsTest(c, 20, mortality2023, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "NewC2023vsA2023.csv")

c <- (mortality2023$SimpleSel == "newA" | mortality2023$SimpleSel == "C") & mortality2023$Sex == "F"
test4 <- statsTest(c, 20, mortality2023, "p-value", newVsOld = TRUE, toLog = TRUE, numDays = 3)
result2 <- extractStats(test4, 3)
write.csv(result2, "NewA2023vsC2023.csv")

