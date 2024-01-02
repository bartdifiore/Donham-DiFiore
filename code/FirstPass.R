#-----------------------------------------------------------------
## Read libraries
#-----------------------------------------------------------------
library(lme4)
library(lmerTest)
library(tidyverse)
library(piecewiseSEM) # This is the workhorse for doing Structural equation modeling. For the full documentation see the online book at <https://jslefche.github.io/sem_book/>. 
library(DHARMa)



#-----------------------------------------------------------------
## Read in data
#-----------------------------------------------------------------

df <- read.csv("data/KelpDataSPULBio_ForSEM.csv") %>% # I'm going to focus on the sheepshead data as biomass for the time being. 
  as_tibble() %>%
  janitor::clean_names() %>% 
  mutate(id = paste(ca_mpa_name_short, year, status, sep = "_")) %>%
  select(id, classcode, density) %>%
  pivot_wider(names_from = classcode, values_from = density) %>%
  separate(id, into = c("site_name", "year", "status"), sep = "_") %>%
  drop_na() %>%
  mutate(status_num = ifelse(status == "MPA", 1, 0))


meta <- read.csv("data/Site_List_All.csv")

l <- length(unique(df$site_name))

x <- unique(df$site_name)

for(i in 1:l) {
  idx <-which(x[i] == df$site_name)
  j <-which(x[i] == meta$CA_MPA_Name_Short)
  m <- length(idx)
  if(meta$MPA_Start[j] == 2005) {
    for (n in 1:m) {
      if (df$year[idx[n]] <=2005) {
        df$time[idx[n]] = 0
      } else if (df$year[idx[n]] == 2006) {
        df$time[idx[n]] = 1
      } else if (df$year[idx[n]] == 2007) {
        df$time[idx[n]] = 2
      } else if (df$year[idx[n]] == 2008) {
        df$time[idx[n]] = 3
      } else if (df$year[idx[n]] == 2009) {
        df$time[idx[n]] = 4
      } else if (df$year[idx[n]] == 2010) {
        df$time[idx[n]] = 5
      } else if (df$year[idx[n]] == 2011) {
        df$time[idx[n]] = 6
      } else if (df$year[idx[n]] == 2012) {
        df$time[idx[n]] = 7
      } else if (df$year[idx[n]] == 2013) {
        df$time[idx[n]] = 8
      } else if (df$year[idx[n]] == 2014) {
        df$time[idx[n]] = 9
      } else if (df$year[idx[n]] == 2015) {
        df$time[idx[n]] = 10
      } else if (df$year[idx[n]] == 2016) {
        df$time[idx[n]] = 11
      } else if (df$year[idx[n]] == 2017) {
        df$time[idx[n]] = 12
      } else if (df$year[idx[n]] == 2018) {
        df$time[idx[n]] = 13
      } else if (df$year[idx[n]] == 2019) {
        df$time[idx[n]] = 14
      } else if (df$year[idx[n]] == 2020) {
        df$time[idx[n]] = 15
      }}
  } else if (meta$MPA_Start[j] == 2012) {
    for (n in 1:m) {
      if (df$year[idx[n]] <= 2012) {
        df$time[idx[n]] = 0
      } else if (df$year[idx[n]] == 2013) {
        df$time[idx[n]] = 1
      } else if (df$year[idx[n]] == 2014) {
        df$time[idx[n]] = 2
      } else if (df$year[idx[n]] == 2015) {
        df$time[idx[n]] = 3
      } else if (df$year[idx[n]] == 2016) {
        df$time[idx[n]] = 4
      } else if (df$year[idx[n]] == 2017) {
        df$time[idx[n]] = 5
      } else if (df$year[idx[n]] == 2018) {
        df$time[idx[n]] = 6
      } else if (df$year[idx[n]] == 2019) {
        df$time[idx[n]] = 7
      } else if (df$year[idx[n]] == 2020) {
        df$time[idx[n]] = 8
      }}}
  else {}
}

#-----------------------------------------------------------------
## DAG 1 (Hypothesis 1): Bottom-up
#-----------------------------------------------------------------

# Kelp -----> purple <-------> red --------> sheephead

bottom_sem <- psem(
  lmer(STRPURAD ~ MACPYRAD + (1|site_name) + (1|year), data = df), 
  lmer(MESFRAAD ~ MACPYRAD + (1|site_name) + (1|year), data = df),
  lmer(SPUL ~ MACPYRAD + STRPURAD + MESFRAAD + (1|site_name) + (1|year), data = df)
)
  
bottom_sem2 <- update(bottom_sem, MESFRAAD %~~% STRPURAD)
summary(bottom_sem2)


# Couldn't get this working yet. Not sure if its implemented for random effects models yet. Seems like it would be possible based on https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.3502 But not sure. Maybe reach out to Jon Lefcheck?  

    # pmultigroup <- piecewiseSEM::multigroup(bottom_sem, group = "status_num")
    # 
    # summary(pmultigroup)
    # pmultigroup
    # 
    # 
    # jutila <- psem(
    #   lm(rich ~ elev + mass, data = meadows),
    #   lm(mass ~ elev, data = meadows)
    # )
    # 
    # jutila.multigroup <- multigroup(jutila, group = "grazed")
    # 
    # jutila.multigroup



#-----------------------------------------------------------------
## DAG 2 (Hypothesis 2): Top-down
#-----------------------------------------------------------------

# Kelp <----- purple <-------> red <-------- sheephead

df.late <- filter(df, time >= 8)


top_sem <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL + status + (1|site_name) + (1|year), data = df.late), 
  lmer(MESFRAAD ~ SPUL + status + (1|site_name) + (1|year), data = df.late),
  lmer(STRPURAD ~ SPUL + status + (1|site_name) + (1|year), data = df.late)
)

top_sem2 <- update(top_sem, MESFRAAD %~~% STRPURAD)
summary(top_sem2, standardize = "scale")


top_sem3 <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL + status + (1|site_name), data = df.late), 
  lmer(MESFRAAD ~ SPUL + status + (1|site_name), data = df.late),
  lmer(STRPURAD ~ SPUL + status + (1|site_name), data = df.late)
)

top_sem4 <- update(top_sem3, MESFRAAD %~~% STRPURAD)
summary(top_sem4, standardize = "scale")

plot(top_sem4)


top_sem0 <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL + (1|site_name) + (1|year), data = df.late), 
  lmer(MESFRAAD ~ SPUL + (1|site_name) + (1|year), data = df.late),
  lmer(STRPURAD ~ SPUL + (1|site_name) + (1|year), data = df.late)
)

top_sem0b <- update(top_sem0, MESFRAAD %~~% STRPURAD)
summary(top_sem0b, standardize = "scale")


pmultigroup <- multigroup(top_sem0b, group = "status")

summary(pmultigroup)
pmultigroup


#-----------------------------------------------------------------
## DAG 3 (Hypothesis 3): Blend
#-----------------------------------------------------------------

# Kelp -----> purple <-------> red <-------- sheephead

blend_sem <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + (1|site_name) + (1|year), data = df), 
  lmer(MESFRAAD ~ SPUL + (1|site_name) + (1|year), data = df),
  lmer(STRPURAD ~ SPUL + (1|site_name) + (1|year), data = df)
)

blend_sem2 <- update(blend_sem, MESFRAAD %~~% STRPURAD)
summary(blend_sem2)


#-----------------------------------------------------------------------------
## Technique based on splitting the data -- Jon Lefcheck does not recommend
#-----------------------------------------------------------------------------


df.mpa <- df %>% filter(status == "mpa")
df.ref <- df %>% filter(status == "reference")

mpa_sem <- psem(
  lmer(MACPYRAD ~ (STRPURAD + MESFRAAD + SPUL)*time + (1|site_name) + (1|year), data = df.mpa), 
  lmer(MESFRAAD ~ (SPUL)*time +  (1|site_name) + (1|year), data = df.mpa),
  lmer(STRPURAD ~ (SPUL)*time + (1|site_name) + (1|year), data = df.mpa)
)

mpa_sem2 <- update(mpa_sem, MESFRAAD %~~% STRPURAD)
summary(mpa_sem2)

ref_sem <- psem(
  lmer(MACPYRAD ~ (STRPURAD + MESFRAAD + SPUL)*time + (1|site_name) + (1|year), data = df.ref), 
  lmer(MESFRAAD ~ (SPUL)*time +  (1|site_name) + (1|year), data = df.ref),
  lmer(STRPURAD ~ (SPUL)*time + (1|site_name) + (1|year), data = df.ref)
)

ref_sem2 <- update(ref_sem, MESFRAAD %~~% STRPURAD)
summary(ref_sem2)

#-----------------------------------------------------------------------------
## Univariate modeling
#-----------------------------------------------------------------------------

mod1 <- lmerTest::lmer(MACPYRAD ~ (status + time)*(scale(STRPURAD) + scale(MESFRAAD) + scale(SPUL)) + (1|site_name) + (1|year), data = df)
summary(mod1)

res <- simulateResiduals(mod1)
plot(res)
hist(residuals(mod1))
anova(mod1)

out1 <- ggeffects::ggpredict(mod1, terms = c("STRPURAD", "status", "time[0,15]" ))
plot(out1)

out2 <- ggeffects::ggpredict(mod1, terms = c("MESFRAAD", "status", "time[0,15]" ))
plot(out2)

out3 <- ggeffects::ggpredict(mod1, terms = c("SPUL", "status", "time[0,15]" ))
plot(out3)


out2 <- ggeffects::ggpredict(mod1, terms = ~SPUL*status)
plot(out2)



out1 <- ggeffects::ggpredict(mod1, terms = c("MESFRAAD", "status"))
plot(out1)














mod2 <- lmer(STRPURAD ~ status*(MESFRAAD + SPUL) + (1|site_name) + (1|year), data = df)
summary(mod2)

out3 <- ggeffects::ggpredict(mod2, terms = ~MESFRAAD*status)
plot(out3)



library(glmmTMB)
mod3 <- glmmTMB(MACPYRAD ~ time + status*(STRPURAD + MESFRAAD + SPUL) + (1|site_name), data = df, family = ziGamma(link = "log"), ziformula = ~1)
summary(mod3)
res <- simulateResiduals(mod3)
plot(res)
