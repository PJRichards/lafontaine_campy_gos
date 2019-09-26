####################################################################
#
# Figure 3. Community alpha-diversity and OTU RA at phyla level
#
####################################################################

# Please note that bar width was standardised using
# code adapted from https://stackoverflow.com/users/8449629/z-lin
# thank you

# packages
library(tidyverse)  # version 1.2.1
library(cowplot)    # version 0.9.4
library(gtable)     # version 0.3.0
library(grid)
library(FSA)        # version 0.8.24


# read in data
tax <- read_tsv("data/mothur/cmpgos.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy") 

OTU <- read_tsv("data/mothur/cmpgos.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared") %>% 
          select(-c(label, numOtus))  %>% 
          filter(str_detect(Group, "kit", negate = TRUE)) %>% 
          mutate(row_sum = rowSums(select(., starts_with("Otu")))) %>%  
          mutate_if(is.numeric, funs(100*(. / row_sum))) %>% 
          select(-row_sum) %>% 
          gather(OTU, value, -Group) %>% 
          spread(Group, value)

alpha <- read_tsv("data/mothur/cmpgos.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.ave-std.summary") %>% 
         filter(group != "e_gos_non_08_1" & group != "e_ctl_cmp_15_5" &
                method == "ave") %>%
         select(-c(method,label), -contains("_lci"), -contains("_hci")) %>% 
         gather(metric, value, -group) %>% 
         mutate(trial = str_sub(group, 1, 1),
                trial = if_else(trial == "e", "early", 
                         if_else(trial == "l", "late", "fail")),
                cohort = str_sub(group, 3, 9),
                cohort = if_else(cohort == "ctl_non", "Control", 
                           if_else(cohort == "ctl_cmp", "Campylobacter",
                             if_else(cohort == "gos_non", "GOS", 
                               if_else(cohort == "gos_cmp", "GOS + Campylobacter", 
                                                  "fail")))),
                cohort = fct_relevel(cohort, "Control", "GOS", "Campylobacter", "GOS + Campylobacter"),
                metric = as.factor(if_else(metric == "invsimpson", "Inverse Simpson", 
                                    if_else(metric == "chao", "Chao",
                                            if_else(metric == "shannon", "Shannon", 
                                                    "fail")))),
                age = str_sub(group, 11, 12) )
            

#################################
#
# alpha-diversity
#
#################################

# alpha-diversity stats

## 6-dc 

### Shannon

#### 8
e8da_shan_anova <- alpha %>%  
                      filter(metric == "Shannon" & 
                              trial == "early" & age == "08") %>% 
                      select(-c(group, metric, trial, age)) %>% 
                      aov(data = ., value ~ cohort) %>% 
                      TukeyHSD(.)  

e8da_shan_anova_label <- as_data_frame(e8da_shan_anova$cohort) %>% 
                            mutate(sig = if_else(.$'p adj' <= 0.001, "***", 
                                                 if_else(.$'p adj' <= 0.01, "**",
                                                         if_else(.$'p adj' <= 0.05, "*", 
                                                                 "ns"))))

#### 15 - 35

e15_35da_shan_ttest <- alpha %>% 
                        mutate(cohort = 
                                 if_else(as.character(cohort) == "GOS + Campylobacter", 
                                         "GOSCampy", as.character(cohort))) %>% 
                        select(-group) %>% 
                        group_by(trial, cohort, age, metric) %>% 
                        summarise(value = list(value)) %>% 
                        group_by(trial, cohort, age, metric) %>% 
                        filter(metric == "Shannon" & 
                                 trial == "early" & age != "08") %>% 
                        spread(cohort, value) %>% 
                        mutate(pvalue = t.test(unlist(Campylobacter), 
                                                unlist(GOSCampy))$p.value,
                               sig = if_else(pvalue <= 0.001, "***", 
                                        if_else(pvalue <= 0.01, "**",
                                          if_else(pvalue <= 0.05, "*", 
                                            "ns"))))

### Inverse Simpson

#### 8
e8da_invsimp_anova <- alpha %>%  
                        filter(metric == "Inverse Simpson" & 
                                 age == "08" & trial == "early") %>% 
                        select(-c(group, metric, trial, age)) %>% 
                        aov(data = ., value ~ cohort) %>% 
                        TukeyHSD(.)   

e8da_invsimp_anova_label <- as_data_frame(e8da_invsimp_anova$cohort) %>% 
                              mutate(sig = if_else(.$'p adj' <= 0.001, "***", 
                                            if_else(.$'p adj' <= 0.01, "**",
                                              if_else(.$'p adj' <= 0.05, "*", 
                                                "ns"))))




#### 15 - 35
# ns
#e15_35da_ttest <- alpha %>% 
#                    mutate(cohort = 
#                              if_else(as.character(cohort) == "GOS + Campylobacter", 
#                                      "GOSCampy", as.character(cohort))) %>% 
#                      select(-group) %>% 
#                      group_by(trial, cohort, age, metric) %>% 
#                      summarise(value = list(value)) %>% 
#                      group_by(trial, cohort, age, metric) %>% 
#                      filter(metric == "Inverse Simpson" & 
#                               trial == "early" &  age != "08") %>% 
#                      spread(cohort, value) %>% 
#                      mutate(pvalue = t.test(unlist(Campylobacter), 
#                                              unlist(GOSCampy))$p.value, 
#                              sig = if_else(pvalue <= 0.001, "***", 
#                                     if_else(pvalue <= 0.01, "**",
#                                      if_else(pvalue <= 0.05, "*", 
#                                        "ns"))))


### Chao

#### 8
# ns
#e8da_chao_kruskal <- alpha %>% 
#                      filter(metric == "Chao" & 
#                             trial == "early" &  age == "08") %>% 
#                      dunnTest(data=., value ~ cohort,  method="bh") 


#e8da_chao_kruskal_label <- as_data_frame(e8da_chao_kruskal$res) %>% 
#                            mutate(sig = if_else(.$'P.adj' <= 0.001, "***", 
#                                          if_else(.$'P.adj' <= 0.01, "**",
#                                            if_else(.$'P.adj' <= 0.05, "*", 
#                                              "ns"))))


#### 15 - 35
## delete "e_gos_cmp_15_2" to match cohort sizes
## delete xxx

#e15_35da_chao_kruskal <- alpha %>% 
#                            filter(group != "e_gos_cmp_15_2" & group != "e_gos_cmp_35_6") %>% 
#                            mutate(cohort = if_else(
#                              as.character(cohort) == "GOS + Campylobacter", 
#                                        "GOSCampy", as.character(cohort))) %>% 
#                            select(-group) %>% 
#                            group_by(trial, cohort, age, metric) %>% 
#                            summarise(value = list(value)) %>% 
#                            group_by(trial, cohort, age, metric) %>% 
#                            filter(metric == "Chao" & 
#                                     trial == "early" &  age != "08") %>% 
#                            spread(cohort, value) %>% 
#                            mutate(pvalue = kruskal.test(unlist(Campylobacter), 
#                                                unlist(GOSCampy))$p.value) %>% 
#                                   sig = if_else(pvalue <= 0.001, "***", 
#                                          if_else(pvalue <= 0.01, "**",
#                                            if_else(pvalue <= 0.05, "*", 
#                                              "ns"))))



## 20-dc 

### Shannon

#### 22
# ns
#l_shan_ttest <- alpha %>% 
#                  mutate(cohort = if_else(as.character(cohort) == "GOS + Campylobacter", 
#                                            "GOSCampy", as.character(cohort))) %>%
#                  select(-group) %>% 
#                  group_by(trial, cohort, age, metric) %>% 
#                  summarise(value = list(value)) %>% 
#                  group_by(trial, cohort, age, metric) %>%
#                  filter(metric == "Shannon" & trial == "late") %>% 
#                  spread(cohort, value) %>% 
#                  mutate(pvalue = t.test(unlist(Campylobacter), 
#                                          unlist(GOSCampy))$p.value,
#                         sig = if_else(pvalue <= 0.001, "***", 
#                                 if_else(pvalue <= 0.01, "**",
#                                   if_else(pvalue <= 0.05, "*", 
#                                       "ns"))))




### Inverse Simpson

#### 22
l_invsimp_ttest <- alpha %>% 
                      mutate(cohort = if_else(as.character(cohort) == "GOS + Campylobacter", 
                                        "GOSCampy", as.character(cohort))) %>%
                      select(-group) %>% 
                      group_by(trial, cohort, age, metric) %>% 
                      summarise(value = list(value)) %>% 
                      group_by(trial, cohort, age, metric) %>%
                      filter(metric == "Inverse Simpson" & trial == "late") %>% 
                      spread(cohort, value) %>% 
                      mutate(pvalue = t.test(unlist(Campylobacter), 
                                              unlist(GOSCampy))$p.value,
                             sig = if_else(pvalue <= 0.001, "***", 
                                      if_else(pvalue <= 0.01, "**",
                                        if_else(pvalue <= 0.05, "*", 
                                          "ns"))))

### Chao

#### 22
# ns
#l_chao_ttest <- alpha %>% 
#                  mutate(cohort = if_else(as.character(cohort) == "GOS + Campylobacter", 
#                                      "GOSCampy", as.character(cohort))) %>%
#                  select(-group) %>% 
#                  group_by(trial, cohort, age, metric) %>% 
#                  summarise(value = list(value)) %>% 
#                  group_by(trial, cohort, age, metric) %>%
#                  filter(metric == "Chao" & trial == "late") %>% 
#                  spread(cohort, value) %>% 
#                  mutate(pvalue = t.test(unlist(Campylobacter), 
#                                          unlist(GOSCampy))$p.value,
#                         sig = if_else(pvalue <= 0.001, "***", 
#                                  if_else(pvalue <= 0.01, "**",
#                                    if_else(pvalue <= 0.05, "*", 
#                                       "ns"))))




## created df to adjust plot spacing to add asterisk
e_alpha_space <- data.frame(cohort = rep("Campylobacter", 2), 
                            metric = c("Inverse Simpson", "Shannon"),
                            value = c(45, 5))

l_alpha_space <- data.frame(cohort = "Campylobacter", 
                            metric = "Inverse Simpson",
                            value = 30)


# plot early alpha-diversity
e_alpha.p <- alpha %>% 
                filter(trial == "early") %>% 
                ggplot(aes(x = cohort, y = value, fill = cohort)) +
                    geom_boxplot(outlier.shape = NA) +
                    geom_point(aes(shape = cohort, fill = cohort), 
                                position = position_jitterdodge(jitter.width = 0.5)) +
                geom_point(data = e_alpha_space, aes(x = cohort, y = value), 
                           colour = "white") +
                scale_fill_manual(values = c("#85E6AE", "#C39BD3", "#999999", "#E69F00")) +
                scale_shape_manual(values = c(22, 23, 21, 25)) +
                facet_grid(cols = vars(age), rows = vars(metric), scale = "free") +
                theme_bw() +
                theme(plot.margin = margin(0.92,0.2,0,0.2, "cm"),
                      legend.position = "none",
                      strip.background = element_rect(fill = "#fff2f4", color = "#000000"),
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      axis.text=element_text(colour="#000000"),
                      axis.text.x = element_text(angle = 45, hjust = 1)) +
                labs(x = "")

# Fix facet so all bars have similar width
# code adapted from https://stackoverflow.com/users/8449629/z-lin

# convert ggplot object to grob object
e_alpha_format.p <- ggplotGrob(e_alpha.p)

# optional: take a look at the grob object's layout
gtable_show_layout(e_alpha_format.p)

# get gtable columns corresponding to the facets
e_alpha_facet.columns <- e_alpha_format.p$layout$l[grepl("panel", 
                                                 e_alpha_format.p$layout$name)]

# get the number of unique x-axis values per facet
e_alpha_x.var <- sapply(ggplot_build(e_alpha.p)$layout$panel_scales_x,
                          function(l) length(l$range$range))


# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
widthz <- e_alpha_format.p$widths[e_alpha_facet.columns] * e_alpha_x.var
widthz$arg1 <- c(2, 2, 4, 2, 2)

e_alpha_format.p$widths[e_alpha_facet.columns] <- widthz

# plot late alpha-diversity
l_alpha_p <- alpha %>% 
                filter(trial == "late") %>% 
                ggplot(aes(x = cohort, y = value, fill = cohort)) +
                       geom_boxplot(outlier.shape = NA) +
                       geom_point(aes(shape = cohort, fill = cohort), 
                                  position = position_jitterdodge(jitter.width = 0.5)) +
                       geom_point(data = l_alpha_space, aes(x = cohort, y = value), 
                                  colour = "white") +
                       scale_fill_manual(values = c("#999999", "#E69F00")) +
                       scale_shape_manual(values = c(21, 25)) +
                       facet_grid(cols = vars(age), rows = vars(metric), scale = "free") +
                       theme_bw() +
                       theme(plot.margin = margin(0.92,3.9,0,0.2, "cm"),
                             legend.position = "none",
                             strip.background = element_rect(fill = "#f2fbf6", color = "#000000"),
                             panel.grid.minor = element_blank(),
                             panel.grid.major = element_blank(),
                             axis.text=element_text(colour="#000000"),
                             axis.text.x = element_text(angle = 45, hjust = 1)) +
                        labs(x = "", y = "")

# add asterisk to plot to indicate significance
e_alpha_sig.p <- ggdraw(e_alpha_format.p) + 
                    draw_label(e15_35da_shan_ttest[3,7][[1]], size=20, x=0.73, y=0.426) +
                    draw_label(e15_35da_shan_ttest[4,7][[1]], size=20, x=0.875, y=0.426) +
                    draw_label(e8da_shan_anova_label$sig[6], size=20, x=0.31, y=0.426) +
                    draw_line(c(0.28, 0.34), 0.42) +
                    draw_label(e8da_invsimp_anova_label$sig[4], size=20, x=0.25, y=0.64) +
                    draw_line(c(0.22, 0.284), 0.634) +
                    draw_label(e8da_invsimp_anova_label$sig[6], size=20, x=0.31, y=0.6) +
                    draw_line(c(0.28, 0.34), 0.594)
  
  

l_alpha_sig.p <- ggdraw(l_alpha_p) +
                    draw_label(l_invsimp_ttest$sig[[1]], size=20, x=0.3, y=0.64)



###################################
#
# Community composition
#
###################################

# format data

## merge OTU and tax datasets
df <- inner_join(tax, OTU, by = "OTU")
    
          
## format as 'tidy'
df_tidy <- df %>% 
            gather("Group", "RA", -c(OTU, Size, Taxonomy)) %>% 
            separate(Taxonomy, 
                    c("kingdom","phylum","class","order","family","genus"), ";") %>% 
            mutate(status = str_sub(Group, 7, 9),
                   trial = str_sub(Group, 1, 1),
                   diet = str_sub(Group, 3, 5),
                   pen = str_sub(Group, str_length(Group)),
                   age = str_sub(Group, 11, 12))
                
            

## summarize at 'Phylum' taxonomic level
df_late_phylum <- df_tidy %>% 
                      filter(trial == 'l') %>% 
                      mutate(phylum = str_sub(phylum, 1, str_length(phylum)-5),
                             diet = if_else(diet == "ctl", "Campylobacter", 
                                      if_else(diet == "gos", "GOS + Campylobacter",                                                                    "fail"))) %>%
                      group_by(diet, pen, phylum) %>% 
                      summarise(abund = sum(RA))
                    

# sanity check
# df_late_phylum %>% group_by(diet, pen) %>% summarise(total = sum(abund)) %>% pull(total) 

df_early_phylum <- df_tidy %>% 
                      filter(trial == 'e') %>%
                      group_by(diet, status, age, pen, phylum) %>% 
                      summarise(abund = sum(RA)) %>% ungroup() %>% 
                      mutate(phylum = str_sub(phylum, 1, str_length(phylum)-5), 
                             cohort = paste0(diet, sep='_', status), 
                             cohort = if_else(cohort == "ctl_non", "Control", 
                                        if_else(cohort == "ctl_cmp", "Campylobacter",
                                          if_else(cohort == "gos_non", "GOS", 
                                            if_else(cohort == "gos_cmp", "GOS + Campylobacter",                                                               "fail")))),
                             cohort = fct_relevel(cohort, "Control", "GOS", 
                                                  "Campylobacter", "GOS + Campylobacter"))

# sanity check
# df_early_phylum %>% group_by(diet, status, age, pen) %>% summarise(total = sum(abund)) %>% pull(total) 


## set phyla (fill) colours
phy_cols <- c("Actinobacteria" = "#c0392b", "Bacteria_unclassified" = "#dc7633", 
              "Firmicutes" = "#52be80", "Proteobacteria" = "#2E86C1")


# plot phylum
early_phy.p <- df_early_phylum %>% 
                    group_by(cohort, age, phylum) %>% 
                    summarise(mean = mean(abund)) %>% 
                    filter(mean > 1) %>% 
                    ggplot(aes(x = cohort, y = mean)) +
                         geom_bar(stat='identity', aes(fill=phylum)) +
                          facet_grid(cols = vars(age), scales = "free_x") +
                         theme_bw() + 
                         theme(strip.background = element_rect(fill = '#fff2f4', color = "#000000"),
                               legend.position = "none",
                               plot.margin = margin(0.92,0.2,0,0.2, "cm"),
                               axis.text.x = element_text(angle = 45, hjust = 1),
                               axis.text=element_text(colour="#000000"),
                               panel.grid.minor = element_blank(),
                               panel.grid.major = element_blank()) +
                               scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
                               scale_fill_manual(values = phy_cols) +
                               labs(x = "", y = "Relative abundance (%)", 
                                    fill = "Phylum")

# Fix facet so all bars have similar width
# code adapted from https://stackoverflow.com/users/8449629/z-lin

# convert ggplot object to grob object
early_phy_format.p <- ggplotGrob(early_phy.p)

# optional: take a look at the grob object's layout
gtable_show_layout(early_phy_format.p)

# get gtable columns corresponding to the facets (5 & 9, in this case)
facet.columns <- early_phy_format.p$layout$l[grepl("panel", 
                                              early_phy_format.p$layout$name)]

# get the number of unique x-axis values per facet (1 & 3, in this case)
x.var <- sapply(ggplot_build(early_phy.p)$layout$panel_scales_x,
                function(l) length(l$range$range))

# change the relative widths of the facet columns based on
# how many unique x-axis values are in each facet
early_phy_format.p$widths[facet.columns] <- early_phy_format.p$widths[facet.columns] * x.var


#grid.draw(early_phy_format.p)



late_phy.p <- df_late_phylum %>% mutate(age = "22") %>% 
                    group_by(diet, phylum, age) %>% 
                    summarise(mean = mean(abund)) %>% 
                    filter(mean > 1) %>% 
                    ggplot(aes(x = diet, y = mean)) +
                          geom_bar(stat = 'identity', aes(fill = phylum)) +
                          facet_wrap(~ age) +
                          theme_bw() +
                          theme(strip.background = element_rect(fill = '#f2fbf6', color = "#000000"),
                                plot.margin = margin(0.92,0.2,0,0.2, "cm"),
                                axis.text.x = element_text(angle = 45, hjust = 1),
                                axis.text = element_text(colour="#000000"),
                                panel.grid.minor = element_blank(),
                                panel.grid.major = element_blank()) +
                                scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
                                scale_fill_manual(values = phy_cols) +
                                labs(x = "", y = "", fill = "Phylum")



# print .tiff
tiff(filename="results/figures/fig_03_diversity.tiff", 
     width = 200, height = 250, units = "mm", res = 300)


plot_grid(e_alpha_sig.p, l_alpha_sig.p,
          early_phy_format.p, late_phy.p,
          rel_widths = c(0.6, 0.38), 
          labels = c("A. 6-dc", "B. 20-dc", "C. 6-dc", "D. 20-dc"), 
          label_x = c(0.042, 0.062, 0.042, 0.062),
          nrow = 2) 

dev.off()


