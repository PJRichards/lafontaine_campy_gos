################################################
#
# Figure S3. LEfSE plots
#
################################################

# Please note that plot levels are arranged in OTU order
# using a neat approach from https://github.com/STAT545-UBC/Discussion/issues/52
# thank you


# packages
library(tidyverse)  # version 1.2.1
library(cowplot)    # version 0.9.4


# read in data

## taxonomy
tax <- read_tsv("data/mothur/cmpgos.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.0.03.cons.taxonomy") %>% 
        separate(Taxonomy, 
                 c("kingdom","phylum","class","order","family","genus"), ";") 

## LEfSE
early_cmp_08_lefse <- "data/mothur/cmpgos.campy.early_08.0.03.pick.0.03.subsample.0.03.lefse_summary"
early_cmp_15_lefse <- "data/mothur/cmpgos.campy.early_15.0.03.pick.0.03.subsample.0.03.lefse_summary"
early_cmp_22_lefse <- "data/mothur/cmpgos.campy.early_22.0.03.pick.0.03.subsample.0.03.lefse_summary"
early_cmp_28_lefse <- "data/mothur/cmpgos.campy.early_28.0.03.pick.0.03.subsample.0.03.lefse_summary"
early_cmp_35_lefse <- "data/mothur/cmpgos.campy.early_35.0.03.pick.0.03.subsample.0.03.lefse_summary"




e_cmp_lefse <- bind_rows(read_tsv(early_cmp_08_lefse),
                          read_tsv(early_cmp_15_lefse),
                            read_tsv(early_cmp_22_lefse),
                              read_tsv(early_cmp_28_lefse),
                                read_tsv(early_cmp_35_lefse)
                                )

e_diet_lefse <- read_tsv("data/mothur/cmpgos.diet.early_08.0.03.pick.0.03.subsample.0.03.lefse_summary")

l_cmp_lefse <- read_tsv("data/mothur/cmpgos.campy.late_22.0.03.pick.0.03.subsample.0.03.lefse_summary")




# format LEfSE
e_cmp_lefse_fmt <- e_cmp_lefse %>% 
                      drop_na() %>%
                      inner_join(tax, by = "OTU") %>% 
                      select(-c(kingdom, phylum, class, order, family)) %>% 
                      mutate(cohort = if_else(str_sub(Class, 3, 5) == "gos", "GOS + Campylobacter", 
                                        if_else(str_sub(Class, 3, 5) == "ctl", "Campylobacter", "fail")),
                             LDA_plot = if_else(cohort == "Campylobacter",  LDA * -1, LDA),
                             ID = paste0(genus %>% str_sub(1, str_length(.)-5), " (", OTU, ")"),
                             comparison = paste0("6-dc (", str_sub(Class, 11, 12)," da)" )) 

l_cmp_lefse_fmt <- l_cmp_lefse %>% 
                      drop_na() %>%
                      inner_join(tax, by = "OTU") %>% 
                      select(-c(kingdom, phylum, class, order, family)) %>% 
                      mutate(cohort = if_else(str_sub(Class, 3, 5) == "gos", "GOS + Campylobacter", 
                                        if_else(str_sub(Class, 3, 5) == "ctl", "Campylobacter", "fail")),
                             LDA_plot = if_else(cohort == "Campylobacter",  LDA * -1, LDA),
                             ID = paste0(genus %>% str_sub(1, str_length(.)-5), " (", OTU, ")"),
                             comparison = paste0("20-dc (", str_sub(Class, 11, 12)," da)" )) 

e_diet_lefse_fmt <- e_diet_lefse %>% 
                      drop_na() %>%
                      inner_join(tax, by = "OTU") %>% 
                      select(-c(kingdom, phylum, class, order, family)) %>% 
                      mutate(cohort = if_else(str_sub(Class, 3, 5) == "ctl", "Campylobacter",
                                        if_else(str_sub(Class, 3, 5) == "gos", "GOS", "fail")),
                             LDA_plot = if_else(cohort == "Campylobacter",  LDA * -1, LDA),
                             ID = paste0(genus %>% str_sub(1, str_length(.)-5), " (", OTU, ")"),
                             comparison = "6-dc (08 da)") 


# plot bars

## describe x-axis labels
lefse_lab = c("-5" = "5", "-4" = "4", "-3" = "3", "-2" = "2", "-1" = "1", 
              "0" = "0", "1" = "1", "2" = "2", "3" = "3", "4" = "4", "5" = "5")

e_cmp_lefse.p <- e_cmp_lefse_fmt %>% 
                        arrange(desc(OTU)) %>%  
                        mutate(ID = factor(ID, unique(ID))) %>% 
                        ggplot(aes(x = ID, y = LDA_plot, fill = cohort)) +
                        geom_bar(stat = "identity", position = "dodge", colour = "black") + 
                        geom_hline(yintercept = 0) +
                        scale_fill_manual(values = c("#999999", "#E69F00")) +
  theme_bw() +
  theme(axis.text = element_text(colour="#000000"),
        axis.text.y = element_text(size = 5),
        strip.background = element_rect(fill = "#fff2f4", color = "#000000"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, colour = "grey50")) +
                        xlab("") + 
                        ylab(expression(paste("LDA score (log")[10]*")")) +
                        scale_y_continuous(breaks=seq(-5, 5, 1), limits=c(-5,5), 
                                                expand=c(0,0), labels=lefse_lab) +
                        coord_flip() +
                        facet_grid(comparison ~ ., scale="free", space="free_y")

l_cmp_lefse.p <- l_cmp_lefse_fmt %>% 
                    arrange(desc(OTU)) %>%  
                    mutate(ID = factor(ID, unique(ID))) %>% 
                    ggplot(aes(x = ID, y = LDA_plot, fill = cohort)) +
                      geom_bar(stat = "identity", position = "dodge", colour = "black") + 
                      geom_hline(yintercept = 0) +
                      scale_fill_manual(values = c("#999999", "#E69F00")) + theme_bw() +
  theme(axis.text = element_text(colour="#000000"),
        axis.text.y = element_text(size = 5),
        strip.background = element_rect(fill = "#f2fbf6", color = "#000000"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2, colour = "grey50")) +
                      xlab("") + 
                      ylab(expression(paste("LDA score (log")[10]*")")) +
                      scale_y_continuous(breaks=seq(-5, 5, 1), limits=c(-5,5), 
                                              expand=c(0,0), labels=lefse_lab) +
                      coord_flip() +
                      facet_grid(comparison ~ .)

e_diet_lefse.p <- e_diet_lefse_fmt %>% 
                      arrange(desc(OTU)) %>%  
                      mutate(ID = factor(ID, unique(ID))) %>% 
                      ggplot(aes(x = ID, y = LDA_plot, fill = cohort)) +
                        geom_bar(stat = "identity", position = "dodge", colour = "black") + 
                        geom_hline(yintercept = 0) +
                        scale_fill_manual(values = c("#999999", "#C39BD3", "#E69F00"),
                                          labels = c("Campylobacter", "GOS", "Campylobacter + GOS"),
                                          limits = c("Campylobacter", "GOS", "Campylobacter + GOS")) +
                        theme_bw() +
                        theme(axis.text = element_text(colour="#000000"),
                              axis.text.y = element_text(size = 5),
                              strip.background = element_rect(fill = "#fff2f4", color = "#000000"),
                              panel.grid.major.x = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.grid.major.y = element_line(size = 0.2, colour = "grey50"),
                              legend.position = "bottom",
                              legend.direction = "vertical") +
                        xlab("") + 
                        ylab(expression(paste("LDA score (log")[10]*")")) +
                        scale_y_continuous(breaks=seq(-5, 5, 1), limits=c(-5,5), 
                                            expand=c(0,0), labels=lefse_lab) +
                        coord_flip() +
                        facet_grid(comparison ~ .)



# draw figure                      
tiff(filename="results/figures/fig_S3_LEfSE.tiff", 
                        width = 225, height = 325, units = "mm", res = 300)

ggdraw() +
        draw_plot(e_cmp_lefse.p + theme(legend.position = "none"), 
                  x = 0, y = 0, width = 0.49, height = 1) +
        draw_plot(l_cmp_lefse.p + theme(legend.position = "none"),
                  x = 0.51, y = 0.89, width = 0.49, height = 0.11) +
        draw_plot(e_diet_lefse.p, 
                  x = 0.51, y = 0.65, width = 0.49, height = 0.21) +
        draw_plot_label(label = c("A", "B", "C"), 
                        size = 15, x = c(0, 0.5, 0.5), y = c(1, 1, 0.86))


dev.off()


