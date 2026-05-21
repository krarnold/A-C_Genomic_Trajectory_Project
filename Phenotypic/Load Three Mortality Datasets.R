# Set data_dir to the path of your local data directory
data_dir <- "."

library(readxl)
library(ggplot2)

mortality <- read.csv(file.path(data_dir, "Mortality_2019.csv"), header = TRUE)
mortality2020 <- read.csv(file.path(data_dir, "Mortality_2020.csv"), header = TRUE)
mortality2023 <- read.csv(file.path(data_dir, "Mortality_2023.csv"), header = TRUE)

numRows <- length(mortality$Alive)
mortality$AlivePerc <- vector(length = numRows)
mortality$AgeSpecificMort <- vector(length = numRows)
mortality$SimpleSel <- vector(length = numRows)
mortality$DeadPerDay <- vector(length = numRows)

mortality$SimpleSel[mortality$Selection == "ACO"] <- "A"
mortality$SimpleSel[mortality$Selection == "AO"] <- "A"
mortality$SimpleSel[mortality$Selection == "CO"] <- "C"
mortality$SimpleSel[mortality$Selection == "NCO"] <- "C"

mortality$SelYear[mortality$SimpleSel == "A"] <- "A2019"
mortality$SelYear[mortality$SimpleSel == "C"] <- "C2019"

mortality$SelPops <- paste(mortality$Selection, mortality$Population)

for (i in c(1:max(mortality$Population))) {
  popAndFemale <- mortality$Population == i & mortality$Sex == "F"
  popAndMale <- mortality$Population == i & mortality$Sex == "M"
  tempMax <- max(mortality$Alive[popAndFemale])
  mortality$AlivePerc[popAndFemale] <- mortality$Alive[popAndFemale] / tempMax
  tempMaxMale <- max(mortality$Alive[popAndMale])
  mortality$AlivePerc[popAndMale] <- mortality$Alive[popAndMale] / tempMaxMale
  
  tempMax <- mortality$Alive[popAndFemale & mortality$Age == 1]
  oldTemp <- tempMax
  newTemp <- oldTemp
  for (j in c(1:max(mortality$Age[mortality$Population == i]))) {
    if (j == 1) {
      mortality$DeadPerDay[popAndFemale & mortality$Age == j] <- 0
      mortality$AgeSpecificMort[popAndFemale & mortality$Age == j] <- 0
    }
    else {
      newTemp <- mortality$Alive[popAndFemale & mortality$Age == j]
      mortality$DeadPerDay[popAndFemale & mortality$Age == j] <- oldTemp - newTemp
      mortality$AgeSpecificMort[popAndFemale & mortality$Age == j] <- (oldTemp - newTemp) / oldTemp
      oldTemp <- newTemp
    }
  }
}

numRows <- length(mortality2020$Alive)
mortality2020$AlivePerc <- vector(length = numRows)
mortality2020$AgeSpecificMort <- vector(length = numRows)
mortality2020$SimpleSel <- vector(length = numRows)
mortality2020$DeadPerDay <- vector(length = numRows)

mortality2020$SimpleSel[mortality2020$Selection == "ACO"] <- "A"
mortality2020$SimpleSel[mortality2020$Selection == "AO"] <- "A"
mortality2020$SimpleSel[mortality2020$Selection == "CO"] <- "C"
mortality2020$SimpleSel[mortality2020$Selection == "NCO"] <- "C"

mortality2020$SelYear[mortality2020$SimpleSel == "A"] <- "A2020"
mortality2020$SelYear[mortality2020$SimpleSel == "C"] <- "C2020"
mortality2020$SelYear[mortality2020$SimpleSel == FALSE] <- FALSE

mortality2020$SimpleSel[mortality2020$Selection == "CACO"] <- "newC"
mortality2020$SimpleSel[mortality2020$Selection == "CAO"] <- "newC"
mortality2020$SimpleSel[mortality2020$Selection == "NACO"] <- "newA"
mortality2020$SimpleSel[mortality2020$Selection == "ANCO"] <- "newA"

mortality2020$SelYear[mortality2020$SimpleSel == "newC"] <- "New C 2020"
mortality2020$SelYear[mortality2020$SimpleSel == "newA"] <- "New A 2020"

mortality2020$SelPops <- paste(mortality2020$Selection, mortality2020$Population)

for (i in c(1:max(mortality2020$Population))) {
  popAndFemale <- mortality2020$Population == i & mortality2020$Sex == "F"
  popAndMale <- mortality2020$Population == i & mortality2020$Sex == "M"
  tempMax <- max(mortality2020$Alive[popAndFemale])
  mortality2020$AlivePerc[popAndFemale] <- mortality2020$Alive[popAndFemale] / tempMax
  tempMaxMale <- max(mortality2020$Alive[popAndMale])
  mortality2020$AlivePerc[popAndMale] <- mortality2020$Alive[popAndMale] / tempMaxMale
  
  tempMax <- mortality2020$Alive[popAndFemale & mortality2020$Age == 1]
  oldTemp <- tempMax
  newTemp <- oldTemp
  for (j in c(1:max(mortality2020$Age[mortality2020$Population == i]))) {
    if (j == 1) {
      mortality2020$DeadPerDay[popAndFemale & mortality2020$Age == j] <- 0
      mortality2020$AgeSpecificMort[popAndFemale & mortality2020$Age == j] <- 0
    }
    else {
      newTemp <- mortality2020$Alive[popAndFemale & mortality2020$Age == j]
      mortality2020$DeadPerDay[popAndFemale & mortality2020$Age == j] <- oldTemp - newTemp
      mortality2020$AgeSpecificMort[popAndFemale & mortality2020$Age == j] <- (oldTemp - newTemp) / oldTemp
      oldTemp <- newTemp
    }
  }
}

numRows <- length(mortality2023$Alive)
mortality2023$AlivePerc <- vector(length = numRows)
mortality2023$AgeSpecificMort <- vector(length = numRows)
mortality2023$SimpleSel <- vector(length = numRows)
mortality2023$DeadPerDay <- vector(length = numRows)

mortality2023$SimpleSel[mortality2023$Selection == "ACO"] <- "A"
mortality2023$SimpleSel[mortality2023$Selection == "AO"] <- "A"
mortality2023$SimpleSel[mortality2023$Selection == "CO"] <- "C"
mortality2023$SimpleSel[mortality2023$Selection == "NCO"] <- "C"

mortality2023$SelYear[mortality2023$SimpleSel == "A"] <- "A2023"
mortality2023$SelYear[mortality2023$SimpleSel == "C"] <- "C2023"
mortality2023$SelYear[mortality2023$SimpleSel == FALSE] <- FALSE

mortality2023$SimpleSel[mortality2023$Selection == "CACO"] <- "newC"
mortality2023$SimpleSel[mortality2023$Selection == "CAO"] <- "newC"
mortality2023$SimpleSel[mortality2023$Selection == "NACO"] <- "newA"
mortality2023$SimpleSel[mortality2023$Selection == "ANCO"] <- "newA"

mortality2023$SelYear[mortality2023$SimpleSel == "newC"] <- "New C 2023"
mortality2023$SelYear[mortality2023$SimpleSel == "newA"] <- "New A 2023"

mortality2023$SelPops <- paste(mortality2023$Selection, mortality2023$Population)

for (i in c(1:max(mortality2023$Population))) {
  popAndFemale <- mortality2023$Population == i & mortality2023$Sex == "F"
  popAndMale <- mortality2023$Population == i & mortality2023$Sex == "M"
  tempMax <- max(mortality2023$Alive[popAndFemale])
  mortality2023$AlivePerc[popAndFemale] <- mortality2023$Alive[popAndFemale] / tempMax
  tempMaxMale <- max(mortality2023$Alive[popAndMale])
  mortality2023$AlivePerc[popAndMale] <- mortality2023$Alive[popAndMale] / tempMaxMale
  
  tempMax <- mortality2023$Alive[popAndFemale & mortality2023$Age == 1]
  oldTemp <- tempMax
  newTemp <- oldTemp
  tempMax2 <- mortality2023$Alive[popAndMale & mortality2023$Age == 1]
  oldTemp2 <- tempMax2
  newTemp2 <- oldTemp2
  for (j in c(1:max(mortality2023$Age[mortality2023$Population == i]))) {
    if (j == 1) {
      mortality2023$DeadPerDay[popAndFemale & mortality2023$Age == j] <- 0
      mortality2023$AgeSpecificMort[popAndFemale & mortality2023$Age == j] <- 0
      mortality2023$DeadPerDay[popAndMale & mortality2023$Age == j] <- 0
      mortality2023$AgeSpecificMort[popAndMale & mortality2023$Age == j] <- 0
    }
    else {
      newTemp <- mortality2023$Alive[popAndFemale & mortality2023$Age == j]
      mortality2023$DeadPerDay[popAndFemale & mortality2023$Age == j] <- oldTemp - newTemp
      mortality2023$AgeSpecificMort[popAndFemale & mortality2023$Age == j] <- (oldTemp - newTemp) / oldTemp
      oldTemp <- newTemp
      newTemp2 <- mortality2023$Alive[popAndMale & mortality2023$Age == j]
      mortality2023$DeadPerDay[popAndMale & mortality2023$Age == j] <- oldTemp2 - newTemp2
      mortality2023$AgeSpecificMort[popAndMale & mortality2023$Age == j] <- (oldTemp2 - newTemp2) / oldTemp2
      oldTemp2 <- newTemp2
    }
  }
}

mortCombined3Years <- rbind(mortality, mortality2020, mortality2023)
