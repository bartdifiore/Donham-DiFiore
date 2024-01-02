library(tidyverse)
library(piecewiseSEM) # This is the workhorse for doing Structural equation modeling. For the full documentation see the online book at <https://jslefche.github.io/sem_book/>. 
library(DHARMa)
library(lme4)
library(nlme)
library(sf)
source("code/theme.R")


#-----------------------------------------------------------------
## Read in data
#-----------------------------------------------------------------

meta <- read.csv("data/Site_List_All.csv") %>% 
  select(Lat, Lon, MPA_Start, CA_MPA_Name_Short)


library(rnaturalearth)
# Just some plotting so I can see it.
      meta_sf <- meta %>% 
        sf::st_as_sf(coords = c("Lon", "Lat"), crs = sf::st_crs(4326))
      
      extent <- st_bbox(meta_sf)
      
      coast <- ne_coastline(scale = 10) %>%
        st_crop(extent)
      
      ggplot()+
        geom_sf(data = meta_sf, color = "darkred")+
        geom_sf(data = coast)

df <- read.csv("data/PISCO.Raw_ish.forSEM.csv") %>%
  as_tibble() %>%
  janitor::clean_names() %>% 
  mutate(ca_mpa_name_short = case_when(ca_mpa_name_short == "Arrow Point to Lion Head Point SMCA" ~ "Arrow Pt to Lion Head Pt SMCA", .default = ca_mpa_name_short)) %>%
  left_join(meta, join_by(ca_mpa_name_short == CA_MPA_Name_Short)) %>%
  mutate(time_since_mpastart = year - MPA_Start )

#--------------------------------------------------------------
## Univariate model
#--------------------------------------------------------------

library(glmmTMB)
mod3 <- glmmTMB(macpyrad ~ site_status*(scale(mesfraad) + scale(strpurad) + scale(spul)) + (1|site) + (1|year), data = df, family = nbinom2(link = "log"), ziformula = ~1)
summary(mod3)
res <- simulateResiduals(mod3)
plot(res)

plot(ggeffects::ggpredict(mod3, terms = ~strpurad*site_status))

temp <- glmmTMB(macpyrad ~ site_status + (1|site) + (1|year), data = df, family = nbinom2(link = "log"), ziformula = ~1)
summary(temp)

pred <- ggeffects::ggpredict(temp, terms = ~site_status)
plot(pred)


mod_df <- df %>% 
  select(-c(spul_bio, spul_legal_bio, allurch)) %>%
  pivot_longer(cols = macpyrad:strpurad)


mod_df %>% group_by(name, site_status) %>% 
  summarize(mean = mean(value))

filter()

mod4 <- glmmTMB(value ~ name*site_status + (1|site) + (1|year), data = mod_df, family = nbinom2(link = "log"), ziformula = ~1)
summary(mod4)
res <- simulateResiduals(mod4)
plot(res)


pred <- as.data.frame(ggeffects::ggpredict(mod4, terms = ~name*site_status))

ggplot(pred, aes(x = group, y = predicted ))+
  geom_linerange(aes(ymin = conf.low, ymax = conf.high, group = group), show.legend = F)+
  geom_point(aes(color = group), size = 2.5, show.legend = F)+
  facet_wrap(~x, scales = "free")+
  theme_bd()

ggsave("figures/mpa_by_trophiclevel.png", device = "png", width = 10, height = 6)


#--------------------------------------------------------------
## SEM
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

temp <- glmmTMB(macpyrad ~ spul_legal + mesfraad + strpurad + panint + site_status + (1|site), data = df, family = nbinom2(link = "log"))

summary(temp)

sem1 <- psem(
glmmTMB(macpyrad ~ spul_legal + mesfraad + strpurad + panint + site_status + (1|site), data = df, family = nbinom2(link = "log"), ziformula = ~site_status), 
  glmmTMB(mesfraad ~ spul_legal + site_status + (1|site) + (1|year), data = df, family = nbinom2(link = "log")),
  glmmTMB(strpurad ~ panint + spul_legal + site_status + (1|site) + (1|year), data = df, family = nbinom2(link = "log")), 
  glmmTMB(spul_legal ~ site_status + (1|site) + (1|year), data = df, family = nbinom2(link = "log")),
  glmmTMB(panint ~ site_status + spul_legal + (1|site) + (1|year), data = df, family = nbinom2(link = "log")),
  mesfraad %~~% strpurad,
  data = df
)
basisSet(sem1)
coefs(sem1)
summary(sem1, conserve = T)


df$status_ordinal <- as.integer(ifelse(df$site_status == "mpa", 1, 0))


df.late <- df %>% 
  filter(case_when(site_status == "mpa" ~ time_since_mpastart > 5, 
                   site_status == "reference" ~ time_since_mpastart > -199))


sem2 <- psem(
  glmmTMB(macpyrad ~ spul_legal + mesfraad + strpurad + panint + status_ordinal + (1|site) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~status_ordinal), 
  glmmTMB(mesfraad ~ spul_legal + status_ordinal + (1|site) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~0),
  glmmTMB(strpurad ~ spul_legal + status_ordinal + (1|site) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~0), 
  glmmTMB(spul_legal ~ status_ordinal + (1|site) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~0),
  glmmTMB(panint ~ status_ordinal + (1|site) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~0),
  mesfraad %~~% strpurad
)
summary(sem2, conserve = T)


temp <- glmmTMB(macpyrad ~ site_status*(scale(spul_legal) + scale(mesfraad) + scale(strpurad) + scale(panint)) + (1|site) + (1|year), data = df, family = nbinom2(link = "log"), ziformula = ~site_status)
summary(temp)

pred <- as.data.frame(ggeffects::ggpredict(temp, terms = ~spul_legal*site_status))
plot(ggeffects::ggpredict(temp, terms = ~strpurad*site_status))

ggplot(pred, aes(x = x, y = predicted ))+
  geom_line(aes(color = group), show.legend = F)+
  theme_bd()

temp <- glmmTMB(macpyrad ~ spul_legal + mesfraad + strpurad + panint + status_ordinal + (1|site) + (1|year), data = df.temp, family = "nbinom1")
family(temp)

temp <- glmer.nb(macpyrad ~ spul_legal + mesfraad + strpurad + panint + status_ordinal + (1|site) + (1|year), data = df.temp)
family(temp)



sem3 <- psem(
  glmer(macpyrad ~  mesfraad + strpurad + spul_legal + panint +  status_ordinal + (1|site/transect) + (1|year), data = df, family = "poisson"), 
  glmer(mesfraad ~ spul_legal + status_ordinal + (1|site/transect) + (1|year), data = df, family = "poisson"),
  glmer(strpurad ~ panint + spul_legal + status_ordinal + (1|site/transect) + (1|year), data = df, family = "poisson"), 
  lmer(spul_legal ~ status_ordinal + (1|site/transect) + (1|year), data = df),
  glmer(panint ~ spul_legal + status_ordinal + (1|site/transect) + (1|year), data = df, family = "poisson"),
  mesfraad %~~% strpurad
)
summary(sem3)





top_sem0_lme <- psem(
  lme(macpyrad ~ spul_legal + mesfraad + strpurad + panint , data = df, random = ~ 1 | site),
  lme(mesfraad ~ spul_legal, data = df, random = ~ 1 | site),
  lme(strpurad ~ panint + spul_legal, data = df, random = ~ 1 | site), 
  mesfraad %~~% strpurad,
  data = df
)
summary(top_sem0_lme)

multigroup(top_sem0_lme, group = "status_ordinal")


df.mpa <- filter(df, site_status == "mpa", time_since_mpastart > 5)
df.ref <- filter(df, site_status == "reference")

kelp.null <- glmer(macpyrad ~ spul_legal + mesfraad + strpurad + panint + (1|site) + (1|year), data = df.mpa, family = "poisson")
ss <- getME(kelp.null, name = c("fixef", "theta"))

sem3.mpa <- psem(
  glmmTMB(macpyrad ~ spul_legal + mesfraad + strpurad + panint + (1|site) + (1|year), data = df.mpa, family = nbinom2(link = "log"), ziformula = ~1), 
  glmmTMB(mesfraad ~ spul_legal + (1|site) + (1|year), data = df.mpa, family = nbinom2(link = "log"), ziformula = ~0),
  glmmTMB(strpurad ~ panint + spul_legal + (1|site) + (1|year), data = df.mpa, family = nbinom2(link = "log"), ziformula = ~0), 
  glmmTMB(panint ~ spul_legal + (1|site) + (1|year), data = df.mpa, family = nbinom2(link = "log"), ziformula = ~0), 
  mesfraad %~~% strpurad
)
summary(sem3.mpa, conserve = T)


sem3.ref <- psem(
  glmmTMB(macpyrad ~ spul_legal + mesfraad + strpurad + panint + (1|site) + (1|year), data = df.ref, family = nbinom2(link = "log"), ziformula = ~1), 
  glmmTMB(mesfraad ~ spul_legal + (1|site) + (1|year), data = df.ref, family = nbinom2(link = "log"), ziformula = ~0),
  glmmTMB(strpurad ~ panint + spul_legal + (1|site) + (1|year), data = df.ref, family = nbinom2(link = "log"), ziformula = ~0), 
  glmmTMB(panint ~ spul_legal + (1|site) + (1|year), data = df.ref, family = nbinom2(link = "log"), ziformula = ~0), 
  mesfraad %~~% strpurad
)
summary(sem3.ref, conserve = T)



ggplot(df, aes(x = panint, y = macpyrad))+
  geom_point(aes(color = site_status))

ggplot(df, aes(x = spul, y = macpyrad))+
  geom_point(aes(color = site_status))

ggplot(df, aes(x = mesfraad, y = macpyrad))+
  geom_point(aes(color = site_status))

ggplot(df, aes(x = strpurad, y = macpyrad))+
  geom_point(aes(color = site_status))







df.late <- df %>% 
  filter(case_when(site_status == "mpa" ~ time_since_mpastart > 5, 
                   site_status == "reference" ~ time_since_mpastart > -199)) %>% 
  mutate(year.fct = factor(year, levels = 1:length(unique(df$year))), 
         site = as.factor(site))


sem.late <- psem(
  glmmTMB(macpyrad ~  mesfraad + strpurad + status_ordinal + (1|site/transect) + (1|year), data = df, family = nbinom2(link = "log"), ziformula = ~.), 
  glmmTMB(mesfraad ~ spul_legal + status_ordinal + (1|site/transect) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~1),
  glmmTMB(strpurad ~ panint + spul_legal + status_ordinal + (1|site/transect) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~1), 
  glmmTMB(spul_legal ~ status_ordinal + (1|site/transect) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~1),
  glmmTMB(panint ~ spul_legal + status_ordinal + (1|site/transect) + (1|year), data = df, family = nbinom1(link = "log"), ziformula = ~1),
  mesfraad %~~% strpurad
)
summary(sem.late, conserve = T)







sem_test <- psem(
  glm(macpyrad ~ spul_legal + mesfraad + strpurad + panint + status_ordinal, data = df, family = "poisson"), 
  glm(mesfraad ~ panint + spul_legal + status_ordinal, data = df, family = "poisson"),
  glm(strpurad ~ spul_legal + status_ordinal, data = df, family = "poisson"), 
  glm(spul_legal ~ status_ordinal, data = df, family = "poisson"),
  glm(panint ~ spul_legal + status_ordinal, data = df, family = "poisson"),
  mesfraad %~~% strpurad
)
summary(sem_test, conserve = T)
