########################################################
#
# Figure S2. Gene expression box plots
#     
########################################################

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
                      mutate(age = as_factor(str_sub(sample_id, 11, 12)),
                            trial = str_sub(sample_id, 1, 1),
                            site = as.factor(if_else(str_sub(target, 1, 1) == "c", "cecum",
                                                     if_else(str_sub(target, 1, 1) == "i", "ileum",
                                                             "fail"))),
                            site = fct_relevel(site, "ileum", "cecum"),
                            target = str_sub(target, 2, str_length(target)-4),
                            target = fct_relevel(target, "IL17A", "IL17F", "IL10", 
                                               "CXCLi1", "CXCLi2", "IL6", "IFNy", "IL1B"),
                            cohort = as_factor(str_sub(sample_id, 3, 9))  ) %>% 
                            filter(cohort != "ctl_moc" &  cohort != "gos_moc")


#### plot ileum data ####

# set figure spacing to accomodate symbols indicating significance
early_spacing <- data.frame(cohort = c("Campylobacter", "Campylobacter", "Campylobacter"),
                           age = c("08", "15", "15"),
                           site = c("ileum", "ileum", "cecum"),
                           target = c("IL17A", "CXCLi2", "IL17F"), 
                           value = c(2.5, 1.5, 1))


late_spacing <- data.frame(cohort = rep("Campylobacter", 4), 
                           age = rep("08", 4),
                           site = c("cecum", "cecum", "ileum", "cecum"),
                           target = c("CXCLi1", "IFNy", "IL6", "IL17F"),
                           value = c(1, 2, 1, 2))



target_names <- c("IL17A" = "IL-17A", "IL17F" = "IL-17F", "IL10" = "IL-10",
                  "CXCLi1" =  "CXCLi1", "CXCLi2" = "CXCLi2", "IL6" = "IL-6", 
                  "IFNy" = paste0("IFN",intToUtf8(947)),  "IL1B" = "IL-1ß")



                     


ge_early.p <- ge_log10dct_tidy %>% 
                      filter(trial == "e") %>% 
                      mutate(cohort = if_else(str_sub(sample_id, 3, 9) == "ctl_cmp", "Campylobacter",
                                        if_else(str_sub(sample_id, 3, 9) == "gos_cmp", "GOS +\nCampylobacter",                                                     "fail"))) %>% 
                      drop_na() %>% 
                      ggplot(aes(x = age, y = value, fill = cohort)) +
                          geom_boxplot(outlier.shape = NA) +
                          geom_point(aes(shape = cohort), size = 1, 
                                     position = position_jitterdodge(jitter.width = 0.2)) +
                          geom_point(data = early_spacing, aes(x = age, y = value), 
                                            fill = "white", colour = "white") +
                          scale_fill_manual(values = c("#999999", "#E69F00")) +
                          scale_shape_manual(values = c(21, 25)) +
                          facet_grid(rows = vars(target), cols = vars(site), scale = "free_y", 
                                     labeller = labeller(target = target_names)) +
                          theme_bw() +
                          theme(strip.background = element_rect(fill = "#fff2f4", color = "#000000"),
                                panel.grid.minor = element_blank(),
                                panel.grid.major = element_blank(),
                                axis.text = element_text(colour="#000000")) +
                          geom_vline(xintercept=c(1.5,2.5,3.5, 4.5), color="#e5e5e5") +
                          ylab(expression(paste("Gene expression ("*log[10],")"))) 


# draw late plot
ge_late.p <- ge_log10dct_tidy %>% 
                  filter(trial == "l") %>% 
                  mutate(cohort = if_else(str_sub(sample_id, 3, 9) == "ctl_cmp", "Campylobacter",
                                      if_else(str_sub(sample_id, 3, 9) == "gos_cmp", "GOS +\nCampylobacter",                                                     "fail"))) %>% 
                  drop_na() %>% 
                  ggplot(aes(x = age, y = value, fill = cohort)) +
                      geom_boxplot(outlier.shape = NA) +
                      geom_point(aes(shape = cohort), size = 1, 
                                 position = position_jitterdodge(jitter.width = 0.2)) +
                      geom_point(data = late_spacing, aes(x = age, y = value), 
                                 fill = "white", colour = "white") +
                      scale_fill_manual(values = c("#999999", "#E69F00")) +
                      scale_shape_manual(values = c(21, 25)) +
                      scale_x_discrete(drop=FALSE) +
                      facet_grid(rows = vars(target), cols = vars(site), scale = "free_y", 
                                 labeller = labeller(target = target_names)) +
                      theme_bw() +
                      theme(strip.background = element_rect(fill = "#f2fbf6", color = "#000000"),
                            panel.grid.minor = element_blank(),
                            panel.grid.major = element_blank(),
                            axis.text = element_text(colour="#000000")) +
                      geom_vline(xintercept=c(1.5,2.5,3.5, 4.5), color="#e5e5e5") +
                      ylab(expression(paste("Gene expression ("*log[10],")"))) 

legend <- get_legend(ge_late.p)

## Perform tests for significance


GE_wilcox <- ge_log10dct_tidy %>% 
                      filter(cohort != "ctl_moc" & cohort != "gos_moc") %>% 
                      group_by(trial, site, cohort, target, age) %>% 
                      summarise(list(value)) %>% rename(values = 'list(value)') %>% 
                      spread(cohort, values) %>% 
                      group_by(trial, site, target, age) %>%
                      mutate(wilcox = wilcox.test(unlist(gos_cmp), unlist(ctl_cmp))$p.value) %>% 
                      group_by(trial, site, age) %>% 
                      mutate(padjust = p.adjust(wilcox, method="BH"),
                             sig = if_else(padjust <= 0.001, "***", 
                                    if_else(padjust <= 0.01, "**",
                                      if_else(padjust <= 0.05, "*", 
                                       "ns"))))

# print .tiff
tiff(filename="results/figures/fig_S2_all_GE.tiff", width = 225, height = 200, units = "mm", res = 300)

plot_grid(ge_early.p + theme(plot.margin = margin(0.92,0.2,0,0.2, "cm"),
                             legend.position = "none"), 
          ge_late.p + theme(plot.margin = margin(0.92,0.2,0,0.2, "cm"),
                            legend.position = "none"), 
          legend, rel_widths = c(0.6, 0.6, 0.2), nrow = 1,
          labels = c("A. 6-dc", "B. 20-dc"), label_x = c(0.07, 0.05)) +
  
          # early ileum CXCLi2 15        
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                    target == "CXCLi2" & age == "15") %>% 
               pull(sig), size = 15, x = 0.118, y = 0.469) +
  
          # early ileum IL-10 15        
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                    target == "IL10" & age == "15") %>% 
               pull(sig), size = 15, x = 0.118, y = 0.69) +
  
          # early ileum IL-10 22        
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                    target == "IL10" & age == "22") %>% 
               pull(sig), size = 15, x = 0.148, y = 0.69) +
          
          # early ileum IL-17A 08        
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                      target == "IL17A" & age == "08") %>% 
                                   pull(sig), size = 15, x = 0.087, y = 0.912) +

          # early ileum IL-17A 15        
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                      target == "IL17A" & age == "15") %>% 
                                   pull(sig), size = 15, x = 0.118, y = 0.912) +
  
          # early ileum IL-17A 22        
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                      target == "IL17A" & age == "22") %>% 
                                   pull(sig), size = 15, x = 0.148, y = 0.912) +
  
          # early ileum IL-17A 28      
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                      target == "IL17A" & age == "28") %>% 
                                   pull(sig), size = 15, x = 0.177, y = 0.912) +
  
          # early ileum IL-17F 22
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                    target == "IL17F" & age == "22") %>% 
               pull(sig), size = 15, x = 0.148, y = 0.801) +
  
          # early ileum IL-17F 28
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "ileum" & 
                                    target == "IL17F" & age == "28") %>% 
               pull(sig), size = 15, x = 0.177, y = 0.801) +
  
          # early cecum CXCLi1 08
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "cecum" & 
                                    target == "CXCLi1" & age == "08") %>% 
               pull(sig), size = 15, x = 0.251, y = 0.578) +
  
           # early cecum IL-10 08
           draw_label(GE_wilcox %>% filter(trial == "e" & site == "cecum" & 
                                    target == "IL10" & age == "08") %>% 
               pull(sig), size = 15, x = 0.251, y = 0.69)  +
  
          # early cecum IL-17A 15        
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "cecum" & 
                                      target == "IL17A" & age == "15") %>% 
                                   pull(sig), size = 15, x = 0.284, y = 0.912) +
  
          # early cecum IL-17A 35
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "cecum" & 
                                      target == "IL17A" & age == "35") %>% 
                                   pull(sig), size = 15, x = 0.376, y = 0.912) +
  
          # early cecum IL-17F 15        
          draw_label(GE_wilcox %>% filter(trial == "e" & site == "cecum" & 
                                    target == "IL17F" & age == "15") %>% 
                                pull(sig), size = 15, x = 0.284, y = 0.801) +
          
          # late ileum IL-6 35
          draw_label(GE_wilcox %>% filter(trial == "l" & site == "ileum" & 
                                      target == "IL6" & age == "35") %>% 
                                   pull(sig), size = 15, x = 0.631, y = 0.358) +
  
          # late cecum CXCLi1 22
          draw_label(GE_wilcox %>% filter(trial == "l" & site == "cecum" & 
                                    target == "CXCLi1" & age == "22") %>% 
                                  pull(sig), size = 15, x = 0.74, y = 0.578) +
  
          # late cecum CXCLi1 28
          draw_label(GE_wilcox %>% filter(trial == "l" & site == "cecum" & 
                                    target == "CXCLi1" & age == "28") %>% 
                                pull(sig), size = 15, x = 0.771, y = 0.578) +
  
          # late cecum IFNy 28
          draw_label(GE_wilcox %>% filter(trial == "l" & site == "cecum" & 
                                   target == "IFNy" & age == "28") %>% 
                                  pull(sig), size = 15, x = 0.771, y = 0.247) +
  
          # late cecum IL-17F 22
          draw_label(GE_wilcox %>% filter(trial == "l" & site == "cecum" & 
                                    target == "IL17F" & age == "22") %>% 
                                  pull(sig), size = 15, x = 0.74, y = 0.801) +
  
          # late cecum IL-17F 28
          draw_label(GE_wilcox %>% filter(trial == "l" & site == "cecum" & 
                                    target == "IL17F" & age == "28") %>% 
                                  pull(sig), size = 15, x = 0.771, y = 0.801) +
  
          # late cecum IL1B 28
          draw_label(GE_wilcox %>% filter(trial == "l" & site == "cecum" & 
                                    target == "IL1B" & age == "28") %>% 
                                pull(sig), size = 15, x = 0.771, y = 0.137)

          
dev.off()


