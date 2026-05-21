library(lme4)
library(car)

pupa2024Melted <- reshape2::melt(pupa2024, measure.vars = c("percDevA", "percDevB", "percDevC"))
c <- pupa2024$SimpleSel == "A" | pupa2024$SimpleSel == "C"
s1_AC <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s1_AC)
Anova(s1_AC)$'Pr(>Chisq)'

c <- pupa2024$SimpleSel == "New A" | pupa2024$SimpleSel == "New C"
s1_AC <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s1_AC)
Anova(s1_AC)$'Pr(>Chisq)'

c <- pupa2024$SimpleSel != "B"
s2_AC <- lmer(value ~ SimpleSel2 * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s2_AC)
Anova(s2_AC)$'Pr(>Chisq)'

c <- pupa2024$Treatment == "NCO" | pupa2024$Treatment == "CO"
s3_AC <- lmer(value ~ Treatment * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'

c <- pupa2024$Treatment == "ACO" | pupa2024$Treatment == "AO"
s3_AC <- lmer(value ~ Treatment * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'

c <- pupa2024$Treatment == "CACO" | pupa2024$Treatment == "CAO"
s3_AC <- lmer(value ~ Treatment * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'

c <- pupa2024$SimpleSel == "New C" | pupa2024$SimpleSel == "C"
s3_AC <- lmer(value ~ OldVsNew * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'

c <- pupa2024$SimpleSel == "New A" | pupa2024$SimpleSel == "A"
s4_AC <- lmer(value ~ OldVsNew * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s4_AC)
Anova(s4_AC)$'Pr(>Chisq)'

c <- pupa2024$Treatment == "ANCO" | pupa2024$Treatment == "NACO"
s3_AC <- lmer(value ~ Treatment * Time + (1 | variable/Replicate), data = pupa2024Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'