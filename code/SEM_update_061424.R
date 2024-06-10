library(tidyverse)
library(piecewiseSEM) # This is the workhorse for doing Structural equation modeling. For the full documentation see the online book at <https://jslefche.github.io/sem_book/>. 
library(DHARMa)
library(lme4)
library(nlme)
library(sf)
source("code/theme.R")



#--------------------------------------------------------------
## DAG
#--------------------------------------------------------------

library(dagitty)
library(ggdag)

coords <- list(
  x = c(Kelp = 0, Purple_urchin = -2, Red_urchin = 2, Sheephead = -1, Status = 0, Lobster = 1),
  y = c(Kelp = 0, Purple_urchin = 3, Red_urchin = 3, Sheephead = 5, Status = 6, Lobster = 5)
)

dag1 <- dagify(Kelp ~ Purple_urchin + Red_urchin + Sheephead + Lobster + Status, 
               Red_urchin ~ Sheephead + Status,
               Purple_urchin ~ Sheephead + Lobster + Status,
               Sheephead ~ Status,
               Lobster ~ Status,
               Red_urchin ~~ Purple_urchin, 
               coords = coords)

ggdag(dag1, 
      text_col = "grey50")+
  theme_dag_blank()



#--------------------------------------------------------------
## Get data
#--------------------------------------------------------------

df <- read.csv("data/PISCO_LTER_KFM_lnRR.csv") %>%
  janitor::clean_names() %>%
  as_tibble()




#--------------------------------------------------------------
## SEM of density
#--------------------------------------------------------------

df1 <- df %>% 
  filter(y %in% c("Panulirus interruptus", "Mesocentrotus franciscanus", "Strongylocentrotus purpuratus", "Macrocystis pyrifera", "Semicossyphus pulcher"), resp == "Den" ) %>%
  select(ca_mpa_name_short, year, y, mpa, reference, time, source) %>%
  pivot_longer(cols = c(mpa, reference), names_to = "status", values_to = "density") %>%
  mutate(density = ifelse(density == 1, 0.9999, density)) %>% # This correction is needed to fit with a zero inflated beta distribution...
  pivot_wider(names_from = y, values_from = density) %>%
  rename(site = ca_mpa_name_short, lob = "Panulirus interruptus", red = "Mesocentrotus franciscanus", purple = "Strongylocentrotus purpuratus", kelp = "Macrocystis pyrifera", sheephead = "Semicossyphus pulcher") %>%
  drop_na() %>% 
  filter(time >= 5)
# names(df1) <- c("site", "year", "time_since", "source", "status", "lob", "red", "purple", "kelp", "sheephead")

summary(df1)

library(glmmTMB)

sem1 <- psem(
  glmmTMB(kelp ~ purple + red + lob + (1|site) + (1|year) + (1|source), ziformula = ~1, data = df1, family = beta_family(link = "logit")), 
  glmmTMB(red ~ lob + sheephead + (1|site) + (1|year) + (1|source), ziformula = ~1, data = df1, family = beta_family(link = "logit")),
  glmmTMB(purple ~ lob + sheephead + (1|site) + (1|year) + (1|source), ziformula = ~1, data = df1, family = beta_family(link = "logit")), 
  # glmmTMB(sheephead ~ 1 + (1|site) + (1|year) + (1|source), ziformula = ~1, data = df1, family = beta_family(link = "logit")),
  # glmmTMB(lob ~ 1 + (1|site) + (1|year) + (1|source), ziformula = ~1, data = df1, family = beta_family(link = "logit")),
  red %~~% purple,
  sheephead %~~% lob
)
basisSet(sem1)
coefs(sem1)
summary(sem1, conserve = T)



sem1 <- psem(
  lmer(kelp ~ purple + red + lob + status + (1|site) + (1|year) + (1|source), data = df1), 
  lmer(red ~ lob + sheephead + status + (1|site) + (1|year) + (1|source), data = df1),
  lmer(purple ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df1), 
  lmer(sheephead ~ status + (1|site) + (1|year) + (1|source), data = df1),
  lmer(lob ~ status + (1|site) + (1|year) + (1|source), data = df1),
  red %~~% purple,
  sheephead %~~% lob
)

summary(sem1, conserve = T)

multigroup()
multigroup(sem1, "status")


df1$status_ordinal <- as.integer(ifelse(df1$status == "mpa", 1, 0))

sem1 <- psem(
  lmer(kelp ~ purple + red + lob + status_ordinal  + (1|site) + (1|year) + (1|source), data = df1), 
  lmer(red ~ sheephead + status_ordinal + (1|site) + (1|year) + (1|source), data = df1),
  lmer(purple ~ lob + sheephead + status_ordinal + (1|site) + (1|year) + (1|source), data = df1), 
  lmer(sheephead ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df1),
  lmer(lob ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df1),
  red %~~% purple
)

summary(sem1)




df.mpa <- df1 %>% filter(status == "mpa")
df.reference <- df1 %>% filter(status == "reference")


sem.mpa <- psem(
  lmer(kelp ~ purple + red + lob + (1|site) + (1|year) + (1|source), data = df.mpa), 
  lmer(red ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df.mpa),
  lmer(purple ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df.mpa), 
  # lmer(sheephead ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df.mpa),
  # lmer(lob ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df.mpa),
  red %~~% purple
  # sheephead %~~% lob
)
summary(sem.mpa)

sem.reference <- psem(
  lmer(kelp ~ purple + red + lob + (1|site) + (1|year) + (1|source), data = df.reference), 
  lmer(red ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df.reference),
  lmer(purple ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df.reference), 
  # lmer(sheephead ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df.mpa),
  # lmer(lob ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df.mpa),
  red %~~% purple
  # sheephead %~~% lob
)
summary(sem.reference)

 

df2 <- df %>% 
  filter(y %in% c("Panulirus interruptus", "Mesocentrotus franciscanus", "Strongylocentrotus purpuratus", "Macrocystis pyrifera", "Semicossyphus pulcher"), resp == "Den" ) %>%
  select(ca_mpa_name_short, year, y, ln_diff, time, source) %>%
  pivot_wider(names_from = y, values_from = ln_diff) %>%
  rename(site = ca_mpa_name_short, lob = "Panulirus interruptus", red = "Mesocentrotus franciscanus", purple = "Strongylocentrotus purpuratus", kelp = "Macrocystis pyrifera", sheephead = "Semicossyphus pulcher") %>%
  drop_na()
# names(df1) <- c("site", "year", "time_since", "source", "status", "lob", "red", "purple", "kelp", "sheephead")

sem1 <- psem(
  lmer(kelp ~ purple + red + lob + (1|site) + (1|year) + (1|source), data = df2), 
  lmer(red ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df2),
  lmer(purple ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df2), 
  # lmer(sheephead ~ 1 + (1|site) + (1|year) + (1|source), data = df1),
  # lmer(lob ~ 1 + (1|site) + (1|year) + (1|source), data = df1),
  red %~~% purple
  # sheephead %~~% lob
)
summary(sem1)







sem1 <- psem(
  lmer(kelp ~ status_ordinal/purple + status_ordinal/red + status_ordinal/lob + status_ordinal/sheephead  + (1|site) + (1|year) + (1|source), data = df1), 
  lmer(red ~ status_ordinal/lob + status_ordinal/sheephead  + (1|site) + (1|year) + (1|source), data = df1),
  lmer(purple ~ status_ordinal/lob + status_ordinal/sheephead + (1|site) + (1|year) + (1|source), data = df1), 
  lmer(sheephead ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df1),
  lmer(lob ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df1),
  red %~~% purple
  # sheephead %~~% lob
)
print.sem1 <- summary(sem1)
print.sem1




sem1 <- psem(
  lmer(kelp ~ status_ordinal*(purple + red + lob + sheephead)  + (1|site) + (1|year) + (1|source), data = df1), 
  lmer(red ~ status_ordinal*(lob + sheephead)  + (1|site) + (1|year) + (1|source), data = df1),
  lmer(purple ~ status_ordinal*(lob + sheephead) + (1|site) + (1|year) + (1|source), data = df1), 
  lmer(sheephead ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df1),
  lmer(lob ~ status_ordinal + (1|site) + (1|year) + (1|source), data = df1),
  red %~~% purple
  # sheephead %~~% lob
)
print.sem1 <- summary(sem1)
print.sem1

fit1 <- lmerTest::lmer(kelp ~ status_ordinal*purple + (1|site) + (1|year) + (1|source), data = df1, REML = F)
summary(fit1)


sem2 <- psem(
  lmer(kelp ~ status_ordinal*purple  + (1|site) + (1|year) + (1|source), data = df1)
)
summary(sem2)



library(ggeffects)
plot(ggpredict(fit1, terms = ~purple*status_numeric))

#--------------------------------------------------------------
## SEM of legal biomass
#--------------------------------------------------------------

df %>% 
  mutate(species_name = case_when(y == "Panulirus interruptus" ~ "lob", 
                                  y == "Mesocentrotus franciscanus" ~ "red", 
                                  y == "Strongylocentrotus purpuratus" ~ "purple",
                                  y == "Macrocystis pyrifera" ~ "kelp",
                                  y == "Semicossyphus pulcher" ~ "sheephead", 
                                  y == "Panulirus interruptus legal" ~ "lob_legal"), 
         temp_id = paste(species_name, resp, sep = "-")) %>% 
  select(ca_mpa_name_short, year, temp_id, mpa, reference, time, source) %>%
  pivot_longer(cols = c(mpa, reference), names_to = "status") %>%
  rename(site = ca_mpa_name_short) %>%
  filter()
  
  filter(y %in% c("Panulirus interruptus", "Mesocentrotus franciscanus", "Strongylocentrotus purpuratus", "Macrocystis pyrifera", "Semicossyphus pulcher"), resp %in% c("Den", "Bio") ) %>%
  select(ca_mpa_name_short, year, y, mpa, reference, time, source) %>%
  pivot_longer(cols = c(mpa, reference), names_to = "status", values_to = "density") %>%
  mutate(density = ifelse(density == 1, 0.9999, density)) %>% # This correction is needed to fit with a zero inflated beta distribution...
  pivot_wider(names_from = y, values_from = density) %>%
  rename(site = ca_mpa_name_short, lob = "Panulirus interruptus", red = "Mesocentrotus franciscanus", purple = "Strongylocentrotus purpuratus", kelp = "Macrocystis pyrifera", sheephead = "Semicossyphus pulcher") %>%
  drop_na() %>% 
  filter(time >= 5)




df %>%
  filter(year == 2017, ca_mpa_name_short == "Campus Point SMCA", time == 5, y == "Panulirus interruptus", resp == "Den")


df_mod <- df %>%
  filter(y %in% c("Panulirus interruptus", "Mesocentrotus franciscanus", "Strongylocentrotus purpuratus", "Macrocystis pyrifera", "Semicossyphus pulcher")) %>% 
  filter(resp == "Den") %>%
  filter(time >= 5) %>%
  mutate(species_name = case_when(y == "Panulirus interruptus" ~ "lob", 
                                  y == "Mesocentrotus franciscanus" ~ "red", 
                                  y == "Strongylocentrotus purpuratus" ~ "purple",
                                  y == "Macrocystis pyrifera" ~ "kelp",
                                  y == "Semicossyphus pulcher" ~ "sheephead")) %>% 
  select(year, time, source, species_name, ln_diff, ca_mpa_name_short) %>% 
  rename(site = ca_mpa_name_short) %>% 
  pivot_wider(names_from = species_name, values_from = ln_diff) %>% 
  drop_na()


sem1 <- psem(
    lmer(kelp ~ purple + red + (1|site) + (1|year) + (1|source), data = df_mod), 
    lmer(red ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df_mod),
    lmer(purple ~ lob + sheephead + (1|site) + (1|year) + (1|source), data = df_mod), 
    # lmer(sheephead ~ 1 + (1|site) + (1|year) + (1|source), data = df_mod),
    # lmer(lob ~ 1 + (1|site) + (1|year) + (1|source), data = df_mod),
    red %~~% purple
  )
summary(sem1)  




fit1 <- lmer(Growth ~ lat*Live + (1|year), data = shipley)
summary(fit1)


Fixed effects:
  Estimate Std. Error t value
(Intercept)  50.8650     6.8442   7.432
lat          -0.1564     0.1094  -1.429
Live         -1.8666     6.9614  -0.268
lat:Live      0.1891     0.1111   1.701

sem2 <- psem(
  lmer(Growth ~ lat*Live + (1|year), data = shipley)
)
summary(sem2)



