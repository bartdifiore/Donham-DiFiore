
lrr <- read.csv("data/PISCOLNRR_ForSEM.csv") %>%
  drop_na(lnSheepDiff)

lrr.late <- lrr %>% filter(time >= 5)



sem_lrr <- psem(
  lmer(lnMacroDiff ~ lnPurpDiff + lnMesDiff + lnSheepDiff + lnLobDiff + (1|CA_MPA_Name_Short) + (1|year), data = lrr.late), 
  lmer(lnMesDiff ~ lnSheepDiff + lnLobDiff + (1|CA_MPA_Name_Short) + (1|year), data = lrr.late),
  lmer(lnPurpDiff ~ lnSheepDiff + lnLobDiff + (1|CA_MPA_Name_Short) + (1|year), data = lrr.late), 
  lnPurpDiff %~~% lnMesDiff,
  data = lrr.late
)



coefs(sem_lrr)


library(lmerTest)
lm1 <- lmer(lnMacroDiff ~ lnPurpDiff + lnMesDiff + lnSheepDiff + lnLobDiff + (1|CA_MPA_Name_Short), data = lrr)
summary(lm1)
anova(lm1)
