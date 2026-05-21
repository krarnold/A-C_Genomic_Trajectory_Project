library(lme4)
library(car)

pupa2025Melted <- reshape2::melt(pupa2025, measure.vars = c("percDevA", "percDevB", "percDevC"))
c <- pupa2025$SimpleSel == "A" | pupa2025$SimpleSel == "C"
s1_AC <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s1_AC)
Anova(s1_AC)$'Pr(>Chisq)'

c <- pupa2025$SimpleSel == "New A" | pupa2025$SimpleSel == "New C"
s1_AC <- lmer(value ~ SimpleSel * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s1_AC)
Anova(s1_AC)$'Pr(>Chisq)'

c <- pupa2025$SimpleSel != "B"
s2_AC <- lmer(value ~ SimpleSel2 * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s2_AC)
Anova(s2_AC)$'Pr(>Chisq)'

c <- pupa2025$Treatment == "NCO" | pupa2025$Treatment == "CO"
s3_AC <- lmer(value ~ Treatment * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'

c <- pupa2025$Treatment == "ACO" | pupa2025$Treatment == "AO"
s3_AC <- lmer(value ~ Treatment * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'

c <- pupa2025$Treatment == "CACO" | pupa2025$Treatment == "CAO"
s3_AC <- lmer(value ~ Treatment * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'

c <- pupa2025$SimpleSel == "New C" | pupa2025$SimpleSel == "C"
s3_AC <- lmer(value ~ OldVsNew * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'

c <- pupa2025$SimpleSel == "New A" | pupa2025$SimpleSel == "A"
s4_AC <- lmer(value ~ OldVsNew * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s4_AC)
Anova(s4_AC)$'Pr(>Chisq)'

c <- pupa2025$Treatment == "ANCO" | pupa2025$Treatment == "NACO"
s3_AC <- lmer(value ~ Treatment * Time + (1 | variable/Replicate), data = pupa2025Melted[c,])
Anova(s3_AC)
Anova(s3_AC)$'Pr(>Chisq)'