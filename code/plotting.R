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
