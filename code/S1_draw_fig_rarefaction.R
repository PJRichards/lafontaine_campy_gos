###############################################
#
# Figure S1. Rarefaction curves
#
# Rarefaction curves indicating sampling effort 
# for each community.
#
##############################################

# packages
library(tidyverse) # version 1.2.1


# list no template controls
drop.ntc <- c("e_kit_neg_1", "e_kit_neg_2", "l_kit_neg")


# read in data
raref <- read_tsv(file="data/mothur/cmpgos.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.rarefaction") %>% 
                      select(numsampled, starts_with("0"))

# format names
names(raref)[2:length(raref)] <- substring(names(raref), 6)[2:length(raref)]

# tidy tibble
# exclude NTC
raref_f <- raref %>% 
                  select(-one_of(drop.ntc)) %>% 
                  gather(key = "sample", value = OTUs, 
                         -numsampled, e_ctl_cmp_08_1:l_gos_cmp_22_6) %>% 
                  drop_na(OTUs)

                                        
raref_f_n <- raref_f %>% 
                mutate(cohort = str_sub(sample, 1, str_length(sample)-2), 
                       cohort = as_factor(str_replace_all(cohort, 
                                 c("e_ctl_cmp_08" = "6-dc. Campylobacter 8 da", 
                                   "e_ctl_cmp_15" = "6-dc. Campylobacter 15 da",
                                   "e_ctl_cmp_22" = "6-dc. Campylobacter 22 da",
                                   "e_ctl_cmp_28" = "6-dc. Campylobacter 28 da",
                                   "e_ctl_cmp_35" = "6-dc. Campylobacter 35 da",
                                   "e_ctl_non_08" = "6-dc. Control 8 da",                  
                                   "e_gos_cmp_08" = "6-dc. GOS + Campylobacter 8 da",
                                   "e_gos_cmp_15" = "6-dc. GOS + Campylobacter 15 da",
                                   "e_gos_cmp_22" = "6-dc. GOS + Campylobacter 22 da",
                                   "e_gos_cmp_28" = "6-dc. GOS + Campylobacter 28 da",
                                   "e_gos_cmp_35" = "6-dc. GOS + Campylobacter 35 da",   
                                   "e_gos_non_08" = "6-dc. GOS 8 da",
                                   "l_ctl_cmp_22" = "20-dc. Campylobacter 22 da",
                                   "l_gos_cmp_22" = "20-dc. GOS + Campylobacter 22 da"))),
                        cohort = fct_relevel(cohort, "6-dc. Campylobacter 8 da", 
                                   "6-dc. Campylobacter 15 da", "6-dc. Campylobacter 22 da", 
                                   "6-dc. Campylobacter 28 da", "6-dc. Campylobacter 35 da", 
                                   "6-dc. GOS + Campylobacter 8 da", "6-dc. GOS + Campylobacter 15 da", 
                                   "6-dc. GOS + Campylobacter 22 da", "6-dc. GOS + Campylobacter 28 da", 
                                   "6-dc. GOS + Campylobacter 35 da", "6-dc. Control 8 da", 
                                   "6-dc. GOS 8 da", "20-dc. Campylobacter 22 da", 
                                   "20-dc. GOS + Campylobacter 22 da"),
                        `bird no.` = str_sub(sample, str_length(sample), str_length(sample)))


# plot data
raref.p <- raref_f_n %>% 
                ggplot() +
                    geom_line(aes(x=numsampled, y = OTUs, group = sample, colour = `bird no.`)) +
                    theme_bw() + 
                    theme(strip.background = element_rect(fill = NA, color = "#000000"),
                          panel.border = element_rect(fill = NA, color = "black"),
                          panel.grid.minor = element_blank()) + 
                    labs(y = "Number of OTUS", x = "Number of reads sampled") +
                    facet_wrap(vars(cohort), ncol = 5)



# print plot
tiff(filename="results/figures/fig_S1_rarefaction.tiff", width=300, height=175, units="mm", res=300)

raref.p

dev.off()




