###############################################################################
#
# Figure 1. Growth performance & Campylobacter colonization
#
###############################################################################

# packages
library(tidyverse)  # version 1.2.1
library(broom)      # version 0.5.2
library(scales)     # version 1.0.0
library(cowplot)    # version 0.9.4


# read mass in data
# and add unique bird IDs
mass <- read_tsv("data/zootechnical/mass.txt") %>% 
          mutate(bird = paste0(group, sep ='_', bird),
                 trial_id = if_else(trial == "e", "A. 6-dc", 
                                    if_else(trial == "l", "B. 20-dc", "fail")),
                 diet = if_else(str_sub(sample_id, 3, 5) == "gos", "GOS + Campylobacter", 
                                if_else(str_sub(sample_id, 3, 5) == "ctl", 
                                        "Campylobacter", "fail"))) 


# ross308 2014 growth performance objectives
# see http://en.aviagen.com/brands/ross/products/ross-308
perf_obj2014 <- tibble(age = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                           16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 
                           29, 30, 31, 32, 33, 34, 35),
                   total_mass = c(42, 57, 73, 91, 111, 134, 160, 189, 220, 256, 
                                  294, 336, 381, 429, 480, 535, 593, 655, 719, 
                                  786, 856, 929, 1004, 1082, 1162, 1244, 1328,
                                  1414, 1501, 1590, 1680, 1771, 1863, 1956, 2050, 
                                  2144),
                   diet = c("target", "target", "target", "target", "target", 
                            "target", "target", "target", "target", "target", 
                            "target", "target", "target", "target", "target", 
                            "target", "target", "target", "target", "target", 
                            "target", "target", "target", "target", "target", 
                            "target", "target", "target", "target", "target", 
                            "target", "target", "target", "target", "target", "target"))

# read in campylobcater count data
raw <- read_tsv("data/zootechnical/counts.txt") %>% 
            drop_na() %>% 
            mutate(age = str_sub(sample_id, 11, 12),
                   trial = str_sub(sample_id, 1, 1),
                   diet = if_else(str_sub(sample_id, 3, 5) == "gos", "GOS + Campylobacter", 
                            if_else(str_sub(sample_id, 3, 5) == "ctl", 
                                    "Campylobacter", "fail")),
                   rawcount = if_else(rawcount == 0, 100, rawcount),
                   logcount = log10(rawcount),
                   trial_id = if_else(trial == "e", "C. 6-dc",
                                if_else(trial == "l", "D. 20-dc", "fail")))


# plot growth rate
grate.p <- mass %>% 
               group_by(trial_id, diet, age) %>% 
               summarise(median = median(total_mass), 
                         sem = sd(total_mass)/sqrt(length(total_mass))) %>% 
               ggplot(aes(x = as.numeric(age), y = median, shape = diet, color = diet)) + 
                  geom_point(aes(fill = diet), colour="#000000") +
                  geom_line(aes(linetype = diet, color = diet), size = 1) +
                  geom_errorbar(aes(ymin = median - sem, ymax = median + sem)) +
                  facet_wrap(~trial_id) +
                  scale_x_continuous(name = "Age (days)", breaks = seq(0, 35, 5)) +
                  scale_y_continuous(name = "Mass (g)", breaks = seq(0, 3000, 1000)) +
                  geom_vline(data=filter(mass, trial_id=="A. 6-dc"), 
                             aes(xintercept=6), linetype="dotted") +                        
                  geom_vline(data=filter(mass,trial_id=="B. 20-dc"), 
                             aes(xintercept=20), linetype="dotted") +
                  geom_line(data = perf_obj2014, aes(x = age, y = total_mass, colour = diet, linetype = diet),
                            size = 1) +
                  theme_bw() + 
                  theme(strip.text.x = element_text(size=14, face = 'bold', hjust = 0),
                        strip.background = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        axis.text=element_text(colour="#000000")) +
                  scale_colour_manual(values = c("#999999", "#E69F00", "#000000")) +
                  scale_fill_manual(values = c("#999999","#E69F00", "#000000")) +
                  scale_shape_manual(values = c(21, 25, NA)) + 
                  scale_linetype_manual(values=c("longdash", "solid", "dotdash"), 
                                        breaks=c("ctl", "gos", "target")) +
                  guides(guide_legend("group"), fill = FALSE, linetype = FALSE)
  
  
  


# plot campylobacter data
count.p <- raw %>% 
              ggplot(aes(x = age, y = logcount, fill = diet)) + 
                    geom_bar(data = raw %>%  filter(age != "08", trial == "e" | age != "22"),
                             position= 'dodge', stat = "summary", fun.y = "mean", color = "black") +
                    geom_point(aes(x = age, shape = diet), 
                                position = position_jitterdodge(jitter.height = 0, 
                                           jitter.width = 0.1, dodge.width=0.9)) + 
                    scale_shape_manual(values = c(21, 25)) +
                    scale_fill_manual(values = c("#999999", "#E69F00")) +
                    
                    ylab(expression(italic("C. jejuni")~textstyle("CFU/g (log"[10]*")"))) +
                    geom_hline(yintercept = 2, linetype = "dashed") +
                    xlab("Age (days)") + 
                    facet_wrap(~trial_id) +
                    theme_bw() + 
                    theme(strip.text.x = element_text(size=14, face = 'bold', hjust = 0),
                          strip.background = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank(),
                          axis.text = element_text(colour="#000000")) +
                    coord_cartesian(ylim = c(2, 10))



# print .tiff image 
tiff(filename="results/figures/fig_1_zootechnical.tiff", width=175, height=200, units="mm", res=300)

plot_grid(grate.p, count.p, nrow = 2)

dev.off()




