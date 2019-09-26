#!/usr/bin/env bash

# download public reads
bash code/get_reads.bash

# plot bird performance & *Campylobacter* loads
code/01_draw_fig_1_zootechnical.R

# plot cytokine expression
code/02_draw_figure_GE_box.R

# download and unpack mothur 1.42.3
bash code/get_mothur_linux_64_1_42_3.bash

# get references for mothur
# generate a customized version of the SILVA reference database that targets the V4 region
# using Silva.seed_v128 and Trainset16_022016
bash code/silva.seed.align.bash

# run mothur through quality control steps
code/mothur/mothur code/get_good_seqs.batch

# calculate sequencing error rate
code/mothur/mothur code/get_error.batch

# get shared OTUs
code/mothur/mothur code/get_shared_otus.batch

# get alpha diversity
code/mothur/mothur "#summary.single(shared=data/mothur/cmpgos.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared, calc=invsimpson-chao-shannon, subsample=19773)"

# plot microbiota data
# alpha diversity
# phyla-level taxonomy
code/03_draw_fig_3_microbiota_plots.R

# get differential OTUs
code/mothur/mothur code/get_differential_OTUs.batch

# get OTU-representative .fasta sequencing
code/mothur/mothur code/get_repfasta.batch

# plot rarefaction curves to determine sampling efficiency
code/S1_draw_fig_rarefaction.R

# plot all GE data
code/S2_draw_fig_all_GE_box.R

# plot differential OTU
code/S3_draw_differential_OTU.R

