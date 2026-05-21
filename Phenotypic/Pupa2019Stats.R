library(lme4)
library(car)

xMelted <- reshape2::melt(x, measure.vars = c("percDevA", "percDevB", "percDevC"))
c <- xMelted$SimpleSel == "A" | xMelted$SimpleSel == "C"
s19_AC <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = xMelted[c,])
Anova(s19_AC)
Anova(s19_AC)$'Pr(>Chisq)'

c <- xMelted$SimpleSel == "New A" | xMelted$SimpleSel == "C"
s19_newAC <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = xMelted[c,])
Anova(s19_newAC)
Anova(s19_newAC)$'Pr(>Chisq)'

c <- xMelted$SimpleSel == "New A" | xMelted$SimpleSel == "A"
s19_newAA <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = xMelted[c,])
Anova(s19_newAA)
Anova(s19_newAA)$'Pr(>Chisq)'

x2Melted <- reshape2::melt(x2, measure.vars = c("percDevA", "percDevB", "percDevC"))
c <- x2Melted$SimpleSel == "A" | x2Melted$SimpleSel == "C"
s19_AC_2 <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = x2Melted[c,])
Anova(s19_AC_2)
Anova(s19_AC_2)$'Pr(>Chisq)'

c <- x2Melted$SimpleSel == "New C" | x2Melted$SimpleSel == "C"
s19_newCC <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = x2Melted[c,])
Anova(s19_newCC)
Anova(s19_newCC)$'Pr(>Chisq)'

c <- x2Melted$SimpleSel == "New C" | x2Melted$SimpleSel == "A"
s19_newCA <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = x2Melted[c,])
Anova(s19_newCA)
Anova(s19_newCA)$'Pr(>Chisq)'
