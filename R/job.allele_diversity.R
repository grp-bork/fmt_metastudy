#/usr/bin/Rscript
################################################################################
#Calculate pairwise SNV distances between all samples for a given species
#
#=> process data for one species at a time
#=> load sample and profile data
#=> load gene content data
#=> load SNV data
#=> calculate pairwise allele distances between all samples
#
#2019-11-14
#sebastian.schmidt@embl.de
################################################################################


################################################################################
################################################################################
# Load Packages
suppressPackageStartupMessages(library("tidyverse", warn.conflicts = F, quietly=T))
library("stringr", warn.conflicts = F, quietly=T)
library("data.table", warn.conflicts = F, quietly=T)
library("Matrix", warn.conflicts = F, quietly = T)

#Make R behave marginally less moronic
options(stringsAsFactors = FALSE)
################################################################################
################################################################################


################################################################################
################################################################################
#Get current species ID
args <- commandArgs(TRUE)
species <- args[2]

#Set folder names
folder.base <- "~/" #=> SET TO APPROPRIATE REPO FOLDER
folder.data <- paste0(folder.base, "data/")
folder.gene_cov <- paste0(folder.data, "pangenomes/per_species.all_FMT/gene.cov/")
folder.gene_snv <- paste0(folder.data, "pangenomes/per_species.all_FMT/snv/", species, "/")

#Point to relevant files
file.data_sample <- paste0(folder.data, "data.sample.RData")
file.data_genome <- paste0(folder.data, "genome_data.combined.Rdata")
file.data_cov.filtered <- paste0(folder.gene_cov, species, ".cov.filtered.RData")
file.data_SNV.filtered <- paste0(folder.gene_snv, species, ".snv_filtered.Rdata")

#Check if the expected files exist, and quit if they don't
if (!file.exists(file.data_cov.filtered)) {writeLines(paste(date(), species, "=> filtered gene coverage file does not exist")); q("no")}
if (!file.exists(file.data_SNV.filtered)) {writeLines(paste(date(), species, "=> SNV file does not exist. Skipping.")); q("no")}

#Load input files
load(file.data_sample)
load(file.data_genome)
load(file.data_cov.filtered)
load(file.data_SNV.filtered)
################################################################################
################################################################################


################################################################################
################################################################################
#Compuate pairwise allele distances
#=> only on positions where both samples have coverage
#=> adjust to total genome size / overlap (obtained via gene lengths)
#=> scale agreement / disagreement by allele at each position
################################################################################
################################################################################
#Sub-select genes that belong to an 'extended core' genome (≥80% prevalence across genomes w/in cluster)
#=> the 80% cutoff is lax to account for incomplete MAGs
if (species %in% c(
  "ref_mOTU_v2_0167",
  "ref_mOTU_v2_0183",
  "ref_mOTU_v2_0199",
  "ref_mOTU_v2_0275",
  "ref_mOTU_v2_0466",
  "ref_mOTU_v2_0561",
  "ref_mOTU_v2_0859",
  "ref_mOTU_v2_0898",
  "ref_mOTU_v2_0899",
  "ref_mOTU_v2_1285",
  "ref_mOTU_v2_1309",
  "ref_mOTU_v2_2805",
  "ref_mOTU_v2_4202",
  "ref_mOTU_v2_4234",
  "ref_mOTU_v2_4268",
  "ref_mOTU_v2_4598",
  "ref_mOTU_v2_4720",
  "ref_mOTU_v2_4875",
  "meta_mOTU_v2_5331",
  "meta_mOTU_v2_5339",
  "meta_mOTU_v2_5437",
  "meta_mOTU_v2_5735",
  "meta_mOTU_v2_5741",
  "meta_mOTU_v2_6063",
  "meta_mOTU_v2_6091",
  "meta_mOTU_v2_6371",
  "meta_mOTU_v2_6834",
  "meta_mOTU_v2_6937",
  "meta_mOTU_v2_7076",
  "meta_mOTU_v2_7082",
  "meta_mOTU_v2_7093",
  "meta_mOTU_v2_7275",
  "meta_mOTU_v2_7449",
  "ANI_AL_95_00218",
  "ANI_AL_95_00220",
  "ANI_AL_95_00223",
  "ANI_AL_95_00229",
  "ANI_AL_95_00230",
  "ANI_AL_95_00231",
  "ANI_AL_95_00232",
  "ANI_AL_95_00241"
)) {
  dat.gene %>% as_tibble(rownames = "gene") %>% arrange(desc(prev.genome)) %>% pull(gene) -> tmp; genes.extended_core <- tmp[1:500]
} else {
  dat.gene %>% as_tibble(rownames = "gene") %>% filter(prev.genome >= 0.8) %>% pull(gene) -> genes.extended_core
}
if (length(genes.extended_core) < 50) {writeLines(paste(date(), species, "=> Less than 50 genes in extended core genome. Skipping.")); q("no")}
size.extended_core <- sum(dat.gene[genes.extended_core, "rep_length"])

#Filter SNV data to relevant genes
dat.snv %>% filter(gene %in% genes.extended_core) -> dat.snv.filtered
if (length(unique(dat.snv.filtered$pos.full)) < 50) {writeLines(paste(date(), species, "=> Less than 50 variant positions in extended core genome. Skipping.")); q("no")}

#Filter coverage data accordingly
dat.ver_cov <- dat.ver_cov[genes.extended_core, ]
dat.hor_cov <- dat.hor_cov[genes.extended_core, ]

#Get total horizontal coverage of extended core genome per sample
dat.hor_cov.per_sample <- colSums(dat.gene[genes.extended_core, "rep_length"] * dat.hor_cov)

#Compute allele diversity at each potential SNV position
dat.snv.filtered %>%
  gather(key = "sample", value = "count", -gene, -pos, -pos.full, -base, -freq.global) %>%
  filter(sample %in% data.sample$sample_alias) %>%
  #Get coverage per position
  group_by(pos.full, sample) %>%
  dplyr::mutate(cov = sum(count)) %>%
  ungroup() %>%
  #Remove unobserved positions (as these are not informative)
  filter(cov > 0) %>%
  #Calculate allele diversity per position
  #=> as InvSimpson of allele frequencies, minus 1 (expected freq)
  dplyr::mutate(
    count_sq = (count / cov) ^ 2
  ) %>%
  group_by(pos.full, sample) %>%
  dplyr::summarise(allele_div = (1 / sum(count_sq)) - 1) %>%
  #Summarise per sampple
  group_by(sample) %>%
  dplyr::summarise(
    n.pos = n(),
    allele_div.total = sum(allele_div)
  ) %>%
  ungroup() %>%
  #Normalise by total n of covered positions
  #=> approximated based on horizontal coverage of extended core genome
  left_join(as_tibble(dat.hor_cov.per_sample, rownames = "sample")) %>%
  mutate(
    hor_cov.total = value,
    allele_div = allele_div.total / value,
    species = species
  ) %>%
  dplyr::select(-value) -> dat.allele_diversity

#Store
save(dat.allele_diversity, file = paste0(folder.base, "results/allele.diversity/allele_diversity.", species, ".RData"))

#Report & quit
writeLines(paste(date(), species, "=> DONE."))
q("no")
################################################################################
################################################################################







###
