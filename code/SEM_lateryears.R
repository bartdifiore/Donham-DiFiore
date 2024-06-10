df <- read.csv("data/PISCO.Raw_ish.forSEM.csv") %>%
  as_tibble() %>%
  janitor::clean_names() %>% 
  mutate(ca_mpa_name_short = case_when(ca_mpa_name_short == "Arrow Point to Lion Head Point SMCA" ~ "Arrow Pt to Lion Head Pt SMCA", .default = ca_mpa_name_short)) %>%
  left_join(meta, join_by(ca_mpa_name_short == CA_MPA_Name_Short)) %>%
  group_by(site, year, site_status, ca_mpa_name_short, MPA_Start) %>%
  summarize(across(spul_bio:strpurad, mean)) %>%
  filter(case_when(MPA_Start == 2005 ~ year >= 2005+5, 
                   MPA_Start == 2012 ~ year >= 2012+5))


#-----------------------------------------------------------

df$status_ordinal <- as.integer(ifelse(df$site_status == "mpa", 1, 0))

sem <- psem(
  glmmTMB(macpyrad ~ spul_legal_bio + mesfraad + strpurad + panint + status_ordinal + (1|site) + (1|year), data = df, family = ziGamma(link = "log"), ziformula = ~status_ordinal), 
  glmmTMB(mesfraad ~ spul_legal_bio + status_ordinal + (1|site) + (1|year), data = df, family = ziGamma(link = "log"), ziformula = ~1),
  glmmTMB(strpurad ~ spul_legal_bio + panint + status_ordinal + (1|site) + (1|year), data = df, family = ziGamma(link = "log"), ziformula = ~1), 
  glmmTMB(spul_legal_bio ~ status_ordinal + (1|site) + (1|year), data = df, family = ziGamma(link = "log"), ziformula = ~1),
  glmmTMB(panint ~ status_ordinal + (1|site) + (1|year), data = df, family = ziGamma(link = "log"), ziformula = ~1),
  mesfraad %~~% strpurad
)
summary(sem, conserve = T)


sem <- psem(
  lmer(macpyrad ~ spul_legal_bio + mesfraad + strpurad + panint + status_ordinal + (1|site) + (1|year), data = df), 
  lmer(mesfraad ~ spul_legal_bio + status_ordinal + (1|site) + (1|year), data = df),
  lmer(strpurad ~ spul_legal_bio + panint + status_ordinal + (1|site) + (1|year), data = df), 
  lmer(spul_legal_bio ~ status_ordinal + (1|site) + (1|year), data = df),
  lmer(panint ~ status_ordinal + (1|site) + (1|year), data = df),
  mesfraad %~~% strpurad
)
summary(sem, direction = c("panint <- mesfraad"))
plot(sem)






sem <- psem(
  glmmTMB(macpyrad ~ spul_legal_bio + mesfraad + strpurad + panint + status_ordinal + (1|site) + (1|year), data = df, family = gaussian(link = "identity")), 
  glmmTMB(mesfraad ~ spul_legal_bio + status_ordinal + (1|site) + (1|year), data = df, family = gaussian(link = "identity")),
  glmmTMB(strpurad ~ spul_legal_bio + panint + status_ordinal + (1|site) + (1|year), data = df, family = gaussian(link = "identity")), 
  glmmTMB(spul_legal_bio ~ status_ordinal + (1|site) + (1|year), data = df, family = gaussian(link = "identity")),
  glmmTMB(panint ~ status_ordinal + (1|site) + (1|year), data = df, family = gaussian(link = "identity")),
  mesfraad %~~% strpurad
)
summary(sem, direction = c("panint <- mesfraad"))

library(DiagrammeR)



coords <- list(
  x = c(Kelp = 0, Red_urchin = 2, Purple_urchin = -2, Sheephead = 1.5, Lobster = -1.5, Status = 0),
  y = c(Kelp = 0,  Red_urchin = 1.5, Purple_urchin = 1.5, Sheephead = 3, Lobster = 3, Status = 4)
)
temp_coords <- as.data.frame(coords)



sem_plot <- plot(sem, return = T)
sem_edf <- get_edge_df(sem_plot) %>%
  filter(id != 13) %>%
  rename(id_external = id)
sem_ndf <- get_node_df(sem_plot) %>% 
  filter(!id %in% c(6, 8)) %>%
  mutate(x = temp_coords$x, 
         y = temp_coords$y) %>%
  mutate(label = case_when(label == "macpyrad" ~ "Kelp", 
                           label == "mesfraad" ~ "Red\nurchin", 
                           label == "strpurad" ~ "Purple\nurchin", 
                           label == "spul_legal_bio" ~ "Sheephead", 
                           label == "panint" ~ "Lobster", 
                           label == "status_ordinal" ~ "Status"), 
         nodes = case_when(nodes == "macpyrad" ~ "Kelp", 
                           nodes == "mesfraad" ~ "Red\nurchin", 
                           nodes == "strpurad" ~ "Purple\nurchin", 
                           nodes == "spul_legal_bio" ~ "Sheephead", 
                           nodes == "panint" ~ "Lobster", 
                           nodes == "status_ordinal" ~ "Status"))




custom_sem <- create_graph() %>%
  add_nodes_from_table(
    table = sem_ndf, label_col = label) %>% 
  add_edges_from_table(
    table = sem_edf,
    from_col = from,
    to_col = to, 
    from_to_map = id_external) %>%
  set_edge_attrs()
  nudge_node_positions_ws()
render_graph(custom_sem)  

custom_sem %>%
  export_graph(
    file_name = "figures/SEM1.png",
    file_type = "PNG", width = 1000, height = 1000
  )


