library(tidyverse)
library(piecewiseSEM) # This is the workhorse for doing Structural equation modeling. For the full documentation see the online book at <https://jslefche.github.io/sem_book/>. 
library(DHARMa)
library(lme4)
library(nlme)



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
  mutate(status_num = ifelse(status == "mpa", 1, 0))


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
## DAG
#-----------------------------------------------------------------

library(dagitty)
library(ggdag)

coords <- list(
  x = c(Kelp = 0, Purple_urchin = -2, Red_urchin = 2, Sheephead = 0, Status = -1.5),
  y = c(Kelp = 0, Purple_urchin = 3, Red_urchin = 3, Sheephead = 5, Status = 5)
)

dag1 <- dagify(Kelp ~ Purple_urchin + Red_urchin + Sheephead + Status, 
               Red_urchin ~ Sheephead + Status,
               Purple_urchin ~ Sheephead + Status,
               Sheephead ~ Status,
               Red_urchin ~~ Purple_urchin, 
               coords = coords)

ggdag(dag1, 
      text_col = "grey50")+
  theme_dag_blank()

#----------------------------------------------------------------
## Fit with piecewiseSEM
#----------------------------------------------------------------

# Kelp <----- purple <-------> red <-------- sheephead

df.late <- filter(df, time >= 5)
df.late$status_num <- as.integer(df.late$status_num)
df$status_num <- as.integer(df$status_num)


# No MPA effect
top_sem0 <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL + (1|site_name) + (1|year), data = df.late), 
  lmer(MESFRAAD ~ SPUL + (1|site_name) + (1|year), data = df.late),
  lmer(STRPURAD ~ SPUL + (1|site_name) + (1|year), data = df.late), 
  MESFRAAD %~~% STRPURAD,
  data = df.late
)
summary(top_sem0)



top_sem0_lme <- psem(
  lme(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL, data = df, random = ~1 | site_name ), 
  lme(MESFRAAD ~ SPUL, data = df, random = ~ 1 | site_name  ),
  lme(STRPURAD ~ SPUL, data = df, random = ~ 1 | site_name ), 
  MESFRAAD %~~% STRPURAD, 
  data = as.data.frame(df.late)
)
summary(top_sem0_lme)

multigroup(top_sem0_lme, group = "status_num")

# MPA effect coded as a categorical variable --- piecewiseSEM used emmeans to estimate marginal mean effects
sem1 <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL + status + (1|site_name) + (1|year), data = df.late), 
  lmer(MESFRAAD ~ SPUL + status + (1|site_name) + (1|year), data = df.late),
  lmer(STRPURAD ~ SPUL + status + (1|site_name) + (1|year), data = df.late), 
  lmer(SPUL ~ status + (1|site_name) + (1|year), data = df.late),
  MESFRAAD %~~% STRPURAD,
  data = df.late
)
coefs(sem1)
summary(sem1)

# This model doesn't include the random effect of year
sem2 <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL + status + (1|site_name), data = df.late), 
  lmer(MESFRAAD ~ SPUL + status + (1|site_name), data = df.late),
  lmer(STRPURAD ~ SPUL + status + (1|site_name), data = df.late), 
  lmer(SPUL ~ status + (1|site_name) + (1|year), data = df.late),
  MESFRAAD %~~% STRPURAD,
  data = df.late
)
summary(sem2)


# This model codes status as a binary where 0 = reference and 1 = MPA. Thus the coefficients for each path represent the expected change in y as x changes from the reference to the MPA.

sem3 <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL + status_num + (1|site_name) + (1|year), data = df.late), 
  lmer(MESFRAAD ~ SPUL + status_num + (1|site_name) + (1|year), data = df.late),
  lmer(STRPURAD ~ SPUL + status_num + (1|site_name) + (1|year), data = df.late), 
  lmer(SPUL ~ status_num + (1|site_name) + (1|year), data = df.late),
  MESFRAAD %~~% STRPURAD,
  data = df.late
)
coefs(sem3)
summary(sem3)

# This model includes an interaction between status and each tropic node

sem4 <- psem(
  lmer(MACPYRAD ~ status_num * (STRPURAD + MESFRAAD + SPUL) + (1|site_name) + (1|year), data = df.late), 
  lmer(MESFRAAD ~ SPUL * status_num + (1|site_name) + (1|year), data = df.late),
  lmer(STRPURAD ~ SPUL * status_num + (1|site_name) + (1|year), data = df.late), 
  lmer(SPUL ~ status_num + (1|site_name) + (1|year), data = df.late),
  MESFRAAD %~~% STRPURAD,
  data = df.late
)
coefs(sem4)


# This includes an interaction between status and year and is based off of all time points
sem5 <- psem(
  lmer(MACPYRAD ~ (time + status_num)*(STRPURAD + MESFRAAD + SPUL) + (1|site_name), data = df), 
  lmer(MESFRAAD ~ (time+ status_num)*SPUL + (1|site_name), data = df),
  lmer(STRPURAD ~ (time+status_num)*SPUL + (1|site_name), data = df), 
  lmer(SPUL ~ status_num + (1|site_name), data = df),
  MESFRAAD %~~% STRPURAD,
  data = df
)
coefs(sem5)


# This includes an interaction with year and is based off of all time points.
sem6 <- psem(
  lmer(MACPYRAD ~ (time)*(STRPURAD + MESFRAAD + SPUL) + (1|site_name), data = df), 
  lmer(MESFRAAD ~ (time)*SPUL + (1|site_name), data = df),
  lmer(STRPURAD ~ (time)*SPUL + (1|site_name), data = df), 
  lmer(SPUL ~ status_num + (1|site_name), data = df),
  MESFRAAD %~~% STRPURAD,
  data = df
)
summary(sem6)



AIC(sem1, sem2, sem3, sem4)



# This is for the early period
df.early <- df %>% filter(time < 2)

sem7 <- psem(
  lmer(MACPYRAD ~ STRPURAD + MESFRAAD + SPUL + status + (1|site_name) + (1|year), data = df.early), 
  lmer(MESFRAAD ~ SPUL + status + (1|site_name) + (1|year), data = df.early),
  lmer(STRPURAD ~ SPUL + status + (1|site_name) + (1|year), data = df.early), 
  lmer(SPUL ~ status + (1|site_name) + (1|year), data = df.early),
  MESFRAAD %~~% STRPURAD,
  data = df.early
)
coefs(sem7)

cbind(coefs(sem1)$Response, coefs(sem1)$Predictor, coefs(sem7)$Estimate, coefs(sem1)$Estimate)
coefs(sem1)



#-----------------------------------------------------------------
## Plot a few of these
#-----------------------------------------------------------------




