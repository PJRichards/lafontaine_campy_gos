###########################################
#
# Figure 2. Gene expression bar charts
#
###########################################

# packages
library(tidyverse)  # version 1.2.1
library(cowplot)    # version 0.9.4
library(FSA)        # version 0.8.24



# read in data
ge_dct <- read_tsv(file = "data/zootechnical/dCT_GE.txt") 


# calculate log10
ge_log10dct <- ge_dct %>% mutate_if(is.numeric, log10)



# tidy format data
# add grouping variables
ge_log10dct_tidy <- ge_log10dct %>% 
  gather(target, value, -sample_id) %>% 
  mutate(age = str_sub(sample_id, 11, 12),
         trial = str_sub(sample_id, 1, 1),
         site = str_sub(target, 1, 1),
         target = str_sub(target, 2, str_length(target)-4),
         target = fct_relevel(target, "IL17A", "IL17F", "IL10", 
                              "CXCLi1", "CXCLi2", "IL6", "IFNy", "IL1B"),
         cohort = as_factor(str_sub(sample_id, 3, 9))  )


#### plot ileum data ####

# set figure spacing to accomodate symbols indicating significance
il_early_spacing <- data.frame(cohort = rep("Campylobacter", 2), 
                               target = c("IL17A", "IL10"),
                               value = c(3, -0))


# draw early ileum plot
ge_early_ileum.p <- ge_log10dct_tidy %>% 
                      filter(site == "i" & age == "08" & 
                             target %in% c("IL17A", "IL17F", "IL10")) %>% 
                      mutate(cohort = if_else(str_sub(sample_id, 3, 9) == "ctl_moc", "Control", 
                                       if_else(str_sub(sample_id, 3, 9) == "ctl_cmp", "Campylobacter",
                                        if_else(str_sub(sample_id, 3, 9) == "gos_moc", "GOS", 
                                          if_else(str_sub(sample_id, 3, 9) == "gos_cmp", "GOS +\nCampylobacter",                                                     "fail")))),
                          cohort = fct_relevel(cohort, "Control", "GOS", 
                              "Campylobacter", "GOS +\nCampylobacter"),
                          id = "6-dc (2 dpi)") %>% 
                      drop_na() %>% 
                      ggplot(aes(x = cohort, y = value, fill = cohort)) +
                        geom_boxplot(outlier.shape = NA) +
                        geom_point(aes(shape = cohort), 
                                   position = position_jitterdodge(jitter.width = 0.6)) +
                        geom_point(data = il_early_spacing, aes(x = cohort, y = value), colour = "white") +
                        scale_fill_manual(values = c("#85E6AE", "#C39BD3", "#999999", "#E69F00")) +
                        scale_shape_manual(values = c(22, 23, 21, 25)) +
                        facet_grid(rows = vars(target), cols = vars(id), scale = "free_y", 
                                    labeller = labeller(target = c("IL17A" = "IL-17A", 
                                                "IL17F" = "IL-17F", "IL10" = "IL-10"))) +
                        theme_bw() +
                        theme(legend.position = "none",
                              strip.background = element_rect(fill = "#fff2f4", color = "#000000"),
                              panel.grid.minor = element_blank(),
                              panel.grid.major = element_blank(),
                              axis.text = element_text(colour="#000000")  )  +
                        ylab(expression(paste("Gene expression ("*log[10],")"))) + xlab("") 


# draw late ileum plot
ge_late_ileum.p <- ge_log10dct_tidy %>% 
  filter(site == "i" & age == "22" & trial == "l",
         target %in% c("IL17A", "IL17F", "IL10")) %>% 
  drop_na() %>% 
  mutate(cohort = if_else(str_sub(sample_id, 3, 9) == "ctl_cmp", "Campylobacter",
                          if_else(str_sub(sample_id, 3, 9) == "gos_cmp", "GOS +\nCampylobacter",                                               "fail")),
         id = "20-dc (2 dpi)") %>% 
  ggplot(aes(x = cohort, y = value, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape = cohort), 
             position = position_jitterdodge(jitter.width = 0.2)) +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  scale_shape_manual(values = c(21, 25)) +
  facet_grid(rows = vars(target), cols = vars(id), scale = "free_y", 
             labeller = labeller(target = c("IL17A" = "IL-17A", 
                                            "IL17F" = "IL-17F", "IL10" = "IL-10"))) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#f2fbf6", color = "#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(colour="#000000"))  +
  ylab(expression(paste("Gene expression ("*log[10],")"))) + xlab("") 


## Perform tests for significance on ileal data ##

# early challenged birds #

# early ileum IL17A
il_early_IL17A_KW <- ge_log10dct_tidy %>% 
  filter(age == "08" & target  == "IL17A", site == "i") %>% 
  dunnTest(value ~ cohort, data=., method="bh") 

# document for plot
il_early_IL17A_ctl_cmp_vs_gos_cmp <- if_else(il_early_IL17A_KW$res[2,4] <= 0.001, "***", 
                                             if_else(il_early_IL17A_KW$res[2,4] <= 0.01, "**",
                                                     if_else(il_early_IL17A_KW$res[2,4] <= 0.05, "*", 
                                                             "ns")))

il_early_IL17A_ctl_moc_vs_gos_cmp <- if_else(il_early_IL17A_KW$res[3,4] <= 0.001, "***", 
                                             if_else(il_early_IL17A_KW$res[3,4] <= 0.01, "**",
                                                     if_else(il_early_IL17A_KW$res[3,4] <= 0.05, "*", 
                                                             "ns")))

il_early_IL17A_ctl_cmp_vs_gos_moc <- if_else(il_early_IL17A_KW$res[4,4] <= 0.001, "***", 
                                             if_else(il_early_IL17A_KW$res[4,4] <= 0.01, "**",
                                                     if_else(il_early_IL17A_KW$res[4,4] <= 0.05, "*", 
                                                             "ns")))

il_early_IL17A_ctl_moc_vs_gos_moc <- if_else(il_early_IL17A_KW$res[5,4] <= 0.001, "***", 
                                             if_else(il_early_IL17A_KW$res[5,4] <= 0.01, "**",
                                                     if_else(il_early_IL17A_KW$res[5,4] <= 0.05, "*", 
                                                             "ns")))

# early ileum IL17F
# no significant differences
#ge_log10dct_tidy %>% 
#                    filter(age == "08" & target  == "IL17F", site == "i") %>% 
#                    dunnTest(value ~ cohort, data=., method="bh")


# early ileum IL10
il_early_IL10_KW <- ge_log10dct_tidy %>% 
  filter(age == "08" & target  == "IL10", site == "i") %>% 
  drop_na() %>% 
  dunnTest(value ~ cohort, data=., method="bh")

# document for plot
il_early_IL10_ctl_moc_vs_gos_cmp <- if_else(il_early_IL10_KW$res[3,4] <= 0.001, "***", 
                                            if_else(il_early_IL10_KW$res[3,4] <= 0.01, "**",
                                                    if_else(il_early_IL10_KW$res[3,4] <= 0.05, "*", 
                                                            "ns")))

il_early_IL10_ctl_moc_vs_gos_moc <- if_else(il_early_IL10_KW$res[5,4] <= 0.001, "***", 
                                            if_else(il_early_IL10_KW$res[5,4] <= 0.01, "**",
                                                    if_else(il_early_IL10_KW$res[5,4] <= 0.05, "*", 
                                                            "ns")))

# late challenged birds #

# all
ileum_late_wilcox <- ge_log10dct_tidy %>% 
  filter(trial == "l" & site == "i") %>% 
  group_by(cohort, target, age) %>% 
  summarise(list(value)) %>% rename(values = 'list(value)') %>% 
  spread(cohort, values) %>% #ungroup() %>% 
  group_by(target, age) %>%
  mutate(wilcox = wilcox.test(unlist(gos_cmp), unlist(ctl_cmp))$p.value) %>% 
  ungroup() %>% 
  mutate(padjust = p.adjust (wilcox, method="BH"),
         sig = if_else(padjust <= 0.001, "***", 
                       if_else(padjust <= 0.01, "**",
                               if_else(padjust <= 0.05, "*", 
                                       "ns"))))



#### plot ceca data ####

# set figure spacing to accomodate symbols indicating significance
ca_early_spacing <- data.frame(cohort = "Campylobacter", 
                               target = "IL10",
                               value = 3)

ca_late_spacing <- data.frame(cohort = "Campylobacter", 
                              target = "IL17F",
                              value = 1)


# draw early cecum plot
ge_early_ceca.p <- ge_log10dct_tidy %>% 
  filter(site == "c" & age == "08" & 
           target %in% c("IL17A", "IL17F", "IL10")) %>%
  drop_na() %>% 
  mutate(cohort = if_else(str_sub(sample_id, 3, 9) == "ctl_moc", "Control", 
                          if_else(str_sub(sample_id, 3, 9) == "ctl_cmp", "Campylobacter",
                                  if_else(str_sub(sample_id, 3, 9) == "gos_moc", "GOS", 
                                          if_else(str_sub(sample_id, 3, 9) == "gos_cmp", "GOS +\nCampylobacter",                                             "fail")))),
         cohort = fct_relevel(cohort, "Control", "GOS", 
                              "Campylobacter", "GOS +\nCampylobacter"),
         id = "6-dc (2 dpi)") %>% 
  ggplot(aes(x = cohort, y = value, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape = cohort), position = position_jitterdodge(jitter.width = 0.6)) +
  geom_point(data = ca_early_spacing, aes(x = cohort, y = value), colour = "white") +
  scale_fill_manual(values = c("#85E6AE", "#C39BD3", "#999999", "#E69F00")) +
  scale_shape_manual(values = c(22, 23, 21, 25)) +
  facet_grid(rows = vars(target), cols = vars(id), scale = "free_y", 
             labeller = labeller(target = c("IL17A" = "IL-17A", 
                                            "IL17F" = "IL-17F", "IL10" = "IL-10"))) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#fff2f4", color = "#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(colour="#000000"))  +
  ylab(expression(paste("Gene expression ("*log[10],")"))) + xlab("")



# draw late cecum plot
ge_late_ceca.p <- ge_log10dct_tidy %>% 
  filter(site == "c" & age == "22" & trial == "l" &
           target %in% c("IL17A", "IL17F", "IL10")) %>% 
  drop_na() %>% 
  mutate(cohort = if_else(str_sub(sample_id, 3, 9) == "ctl_cmp", "Campylobacter",
                          if_else(str_sub(sample_id, 3, 9) == "gos_cmp", "GOS +\nCampylobacter",                                               "fail")),
         id = "20-dc (2 dpi)") %>% 
  ggplot(aes(x = cohort, y = value, fill = cohort)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(shape = cohort), 
             position = position_jitterdodge(jitter.width = 0.2)) +
  geom_point(data = ca_late_spacing, aes(x = cohort, y = value), colour = "white") +
  scale_fill_manual(values = c("#999999", "#E69F00")) +
  scale_shape_manual(values = c(21, 25)) +
  facet_grid(rows = vars(target), cols = vars(id), scale = "free_y", 
             labeller = labeller(target = c("IL17A" = "IL-17A", 
                                            "IL17F" = "IL-17F", "IL10" = "IL-10"))) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "#f2fbf6", color = "#000000"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(colour="#000000"))  +
  ylab(expression(paste("Gene expression ("*log[10],")"))) + xlab("") 







## Perform tests for significance on cecal data ##

# early challenged birds #

# early cecum IL17A
# no significant differences
# ge_log10dct_tidy %>% 
#                 filter(age == "08" & target  == "IL17A", site == "c") %>% 
#                 dunnTest(value ~ cohort, data=., method="bh")


# early cecum IL17F
# no significant differences
# ge_log10dct_tidy %>% 
#                 filter(age == "08" & target  == "IL17F", site == "c") %>% 
#                 dunnTest(value ~ cohort, data=., method="bh")


# early cecum IL10
ca_early_IL10_KW <- ge_log10dct_tidy %>% 
  filter(age == "08" & target  == "IL10", site == "c") %>% 
  drop_na() %>% 
  dunnTest(value ~ cohort, data=., method="bh")

# document for plot
ca_early_IL10_ctl_cmp_vs_gos_cmp <- if_else(ca_early_IL10_KW$res[2,4] <= 0.001, "***", 
                                            if_else(ca_early_IL10_KW$res[2,4] <= 0.01, "**",
                                                    if_else(ca_early_IL10_KW$res[2,4] <= 0.05, "*", 
                                                            "ns")))

ca_early_IL10_ctl_moc_vs_gos_cmp <- if_else(ca_early_IL10_KW$res[3,4] <= 0.001, "***", 
                                            if_else(ca_early_IL10_KW$res[3,4] <= 0.01, "**",
                                                    if_else(ca_early_IL10_KW$res[3,4] <= 0.05, "*", 
                                                            "ns")))

ca_early_IL10_ctl_cmp_vs_gos_moc <- if_else(ca_early_IL10_KW$res[4,4] <= 0.001, "***", 
                                            if_else(ca_early_IL10_KW$res[4,4] <= 0.01, "**",
                                                    if_else(ca_early_IL10_KW$res[4,4] <= 0.05, "*", 
                                                            "ns")))

ca_early_IL10_ctl_moc_vs_gos_moc <- if_else(ca_early_IL10_KW$res[5,4] <= 0.001, "***", 
                                            if_else(ca_early_IL10_KW$res[5,4] <= 0.01, "**",
                                                    if_else(ca_early_IL10_KW$res[5,4] <= 0.05, "*", 
                                                            "ns")))



# late challenged bird
cecum_late_wilcox <- ge_log10dct_tidy %>% 
  filter(trial == "l" & site == "c") %>% 
  group_by(cohort, target, age) %>% 
  summarise(list(value)) %>% 
  rename(values = 'list(value)') %>% 
  spread(cohort, values) %>% 
  group_by(target, age) %>%
  mutate(wilcox = wilcox.test(unlist(gos_cmp), 
                              unlist(ctl_cmp))$p.value) %>% 
  ungroup() %>% 
  mutate(padjust = p.adjust (wilcox, method="BH"),
         sig = if_else(padjust <= 0.001, "***", 
                       if_else(padjust <= 0.01, "**",
                               if_else(padjust <= 0.05, "*", 
                                       "ns"))))

# document for plot
ca_late_IL17F_ctl_cmp_vs_gos_cmp <- cecum_late_wilcox[4,7][[1]]


# draw site co-plots
ileum.p <- plot_grid(ge_early_ileum.p + theme(plot.margin = margin(0.92,0.2,0,0.2, "cm")), 
                     ge_late_ileum.p + theme(plot.margin = margin(0.92,0.2,0,0.2, "cm")),
                     rel_widths = c(0.7 , 0.38), 
                     labels = "A. Ileum", label_x = c(0.03, 0.095),
                     nrow = 1)

ceca.p <- plot_grid(ge_early_ceca.p + theme(plot.margin = margin(0.92,0.2,0,0.2, "cm")), 
                    ge_late_ceca.p + theme(plot.margin = margin(0.92,0.2,0,0.2, "cm")),
                    rel_widths = c(0.7 , 0.38), 
                    labels = "B. Cecum", label_x = c(0.03, 0.095),
                    nrow = 1)



# print .tiff
tiff(filename="results/figures/fig_2_GE.tiff", width = 210, height = 210, units = "mm", res = 300)

plot_grid(ileum.p, ceca.p, nrow = 2) +
  ### ileum early
  
  ## IL-17A
  # ctl_cmp vs gos_cmp
  draw_line(x = c(0.4, 0.535), y = 0.896, size = 0.5) +
  draw_label(il_early_IL17A_ctl_cmp_vs_gos_cmp, size = 15, x = 0.47, y = 0.899) +
  
  # ctl_moc_vs_gos_cmp
  draw_line(x = c(0.138, 0.535), y = 0.912, size = 0.5) +
  draw_label(il_early_IL17A_ctl_moc_vs_gos_cmp, size=15, x=0.35, y=0.915) +
  
  # ctl_moc_vs_gos_moc
  draw_line(x = c(0.138, 0.27), y = 0.884, size = 0.5) +
  draw_label(il_early_IL10_ctl_moc_vs_gos_moc, size = 15, x = 0.21, y = 0.887) +
  
  # ctl_cmp_vs_gos_moc
  draw_line(x = c(0.27, 0.4), y = 0.89, size = 0.5) +
  draw_label(il_early_IL17A_ctl_cmp_vs_gos_moc, size = 15, x = 0.342, y = 0.893) +
  
  ## IL-10
  # ctl_moc_vs_gos_cmp 
  draw_line(x = c(0.138, 0.535), y = 0.662, size = 0.5) +
  draw_label(il_early_IL10_ctl_moc_vs_gos_cmp , size=15, x = 0.35, y = 0.665) +
  
  # ctl_moc_vs_gos_moc 
  draw_line(x = c(0.138, 0.27), y = 0.646, size = 0.5) +
  draw_label(il_early_IL10_ctl_moc_vs_gos_moc, size=15, x = 0.21, y = 0.649) +
  
  ### cecum early
  
  ## IL-10
  # ctl_cmp vs gos_cmp
  draw_line(x = c(0.4, 0.535), y = 0.125, size = 0.5) +
  draw_label(ca_early_IL10_ctl_cmp_vs_gos_cmp , size = 15, x = 0.47, y = 0.128) +
  
  # ctl_moc_vs_gos_cmp
  draw_line(x = c(0.138, 0.535), y = 0.162, size = 0.5) +
  draw_label(ca_early_IL10_ctl_moc_vs_gos_cmp, size = 15, x = 0.35, y = 0.165) +
  
  # ctl_moc_vs_gos_moc 
  draw_line(x = c(0.138, 0.27), y = 0.146, size = 0.5) +
  draw_label(ca_early_IL10_ctl_moc_vs_gos_moc , size=15, x = 0.21, y = 0.149) +
  
  # ctl_cmp_vs_gos_moc 
  draw_line(x = c(0.27, 0.4), y = 0.14, size = 0.5) +
  draw_label(ca_early_IL10_ctl_cmp_vs_gos_moc, size = 15, x = 0.342, y = 0.143) +
  
  ### cecum late 
  
  ## IL-17F
  # ctl_cmp vs gos_cmp
  draw_line(x = c(0.785, 0.895), y = 0.285, size = 0.5) +
  draw_label(ca_late_IL17F_ctl_cmp_vs_gos_cmp, size = 15, x = 0.84, y = 0.288) 


dev.off()


