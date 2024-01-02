




sem_saturated <- list()
sem_saturated[[1]] <-  glmmTMB(macpyrad ~ site_status/spul_legal + site_status/mesfraad + site_status/strpurad + site_status/panint + (1|site) + (1|year), data = df, family = nbinom2(link = "log"), ziformula = ~site_status) 
sem_saturated[[2]] <- glmmTMB(mesfraad ~ site_status/spul_legal  + site_status/panint + site_status/strpurad + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))
sem_saturated[[3]] <- glmmTMB(strpurad ~ site_status/spul_legal  + site_status/panint + site_status/mesfraad + (1|site) + (1|year), data = df, family = nbinom2(link = "log")) 
sem_saturated[[4]] <- glmmTMB(spul_legal ~ site_status/panint + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))
sem_saturated[[5]] <- glmmTMB(panint ~ site_status/spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))

get_df <- function(x){
  extractAIC(x)[1]
}
logLik.sat <- sum(unlist(lapply(sem_saturated, logLik)))
df.sat <- sum(unlist(lapply(sem_saturated, get_df)))


get_df(sem_saturated[[1]])

sem_hypothesis <- list()
sem_hypothesis[[1]] <-  glmmTMB(macpyrad ~ site_status/spul_legal + site_status/mesfraad + site_status/strpurad + site_status/panint + (1|site) + (1|year), data = df, family = nbinom2(link = "log"), ziformula = ~site_status) 
sem_hypothesis[[2]] <- glmmTMB(mesfraad ~ site_status/spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))
sem_hypothesis[[3]] <- glmmTMB(strpurad ~ site_status/panint + site_status/spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log")) 
sem_hypothesis[[4]] <- glmmTMB(spul_legal ~ site_status + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))
sem_hypothesis[[5]] <- glmmTMB(panint ~ site_status/spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))

logLik.hypothesis <- sum(unlist(lapply(sem_hypothesis, logLik)))
df.hypothesis <- sum(unlist(lapply(sem_hypothesis, get_df)))


sem_nostatus <- list()
sem_nostatus[[1]] <-  glmmTMB(macpyrad ~ spul_legal + mesfraad + strpurad + panint + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))
sem_nostatus[[2]] <- glmmTMB(mesfraad ~ spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))
sem_nostatus[[3]] <- glmmTMB(strpurad ~ panint + spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log")) 
sem_nostatus[[4]] <- glmmTMB(panint ~ spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))

logLik.nostatus <- sum(unlist(lapply(sem_nostatus, logLik)))
df.nostatus <- sum(unlist(lapply(sem_nostatus, get_df)))



sem_nointeraction <- list()
sem_nointeraction[[1]] <-  glmmTMB(macpyrad ~ site_status + spul_legal + mesfraad + strpurad + panint + (1|site) + (1|year), data = df, family = nbinom2(link = "log"), ziformula = ~site_status)
sem_nointeraction[[2]] <- glmmTMB(mesfraad ~ site_status + spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))
sem_nointeraction[[3]] <- glmmTMB(strpurad ~ site_status + panint + spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log")) 
sem_nointeraction[[4]] <- glmmTMB(spul_legal ~ site_status + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))
sem_nointeraction[[5]] <- glmmTMB(panint ~ site_status + spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log"))

logLik.nointeraction <- sum(unlist(lapply(sem_nointeraction, logLik)))
df.nointeraction <- sum(unlist(lapply(sem_nointeraction, get_df)))


# Model comparison

  #Saturated model agaist the hypothsis
  Chi.hypothsis.saturated <- -2*(logLik.hypothesis - logLik.sat)
  1-pchisq(Chi.hypothsis.saturated,(df.sat-df.hypothesis))

  
  #Saturated model against the nointeraction
  Chi.nointeraction.saturated <- -2*(logLik.nointeraction - logLik.sat)
  1-pchisq(Chi.nointeraction.saturated,(df.sat - df.nointeraction))
  







