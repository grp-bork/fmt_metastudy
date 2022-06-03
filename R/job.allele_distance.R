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
  dat.gene %>% as_tibble(rownames = "gene") %>% arrange(desc(prev.genome)) %>% pull(gene) -> tmp; genes.extended_core <- tmp[1:100]
} else {
  dat.gene %>% as_tibble(rownames = "gene") %>% filter(prev.genome >= 0.8) %>% pull(gene) -> genes.extended_core
}
if (length(genes.extended_core) < 50) {writeLines(paste(date(), species, "=> Less than 50 genes in extended core genome. Skipping.")); q("no")}
size.extended_core <- sum(dat.gene[genes.extended_core, ]$rep_length)

#Filter SNV data to relevant genes
dat.snv %>% filter(gene %in% genes.extended_core) -> dat.snv.filtered
if (length(unique(dat.snv.filtered$pos.full)) < 50) {writeLines(paste(date(), species, "=> Less than 50 variant positions in extended core genome. Skipping.")); q("no")}

#Filter coverage data accordingly
dat.ver_cov <- dat.ver_cov[genes.extended_core, ]
dat.hor_cov <- dat.hor_cov[genes.extended_core, ]

#Precalculate per-sample coverage at each potential SNV position
dat.snv.filtered %>%
  dplyr::select(-freq.global) %>%
  gather(key = "sample", value = "count", -gene, -pos, -pos.full, -base) %>%
  filter(sample %in% data.sample$sample_alias) %>%
  group_by(pos.full, sample) %>%
  dplyr::mutate(cov = sum(count)) %>%
  ungroup() %>%
  filter(cov > 0) %>%
  dplyr::mutate(count = count / cov) %>%
  dplyr::mutate(allele = str_c(pos.full, base, sep=":")) %>%
  dplyr::mutate(idx.sample = match(sample, data.sample$sample_alias)) -> dat.allele
#Add position index
pos.retained <- dat.allele %>% dplyr::select(gene, pos, pos.full) %>% distinct() %>% arrange(gene, pos) %>% pull(pos.full) %>% unique()
allele.retained <- dat.allele %>% dplyr::select(gene, pos, allele) %>% distinct() %>% arrange(gene, pos) %>% pull(allele) %>% unique()
dat.allele %>% mutate(idx.pos.full = match(pos.full, pos.retained), idx.allele = match(allele, allele.retained)) -> dat.allele

#Allocate matrices for alleles and coverage
mat.allele <- sparseMatrix(
  i = dat.allele$idx.allele,
  j = dat.allele$idx.sample,
  x = dat.allele$count,
  dims = c(length(allele.retained), nrow(data.sample)),
  dimnames = list(allele.retained, data.sample$sample_alias)
)
mat.cov <- sparseMatrix(
  i = dat.allele$idx.allele,
  j = dat.allele$idx.sample,
  x = dat.allele$idx.pos.full,
  dims = c(length(allele.retained), nrow(data.sample)),
  dimnames = list(allele.retained, data.sample$sample_alias)
)

#Preallocate results matrices
mat.diff <- mat.shared_cov <- matrix(NA, nrow = nrow(data.sample), ncol = nrow(data.sample), dimnames = list(data.sample$sample_alias, data.sample$sample_alias))

#Iterate through samples and compute pairwise allele distances
for (s in data.sample$sample_alias) {
  #Subset data to positions observed in current focal sample
  keep.alleles <- mat.cov[, s] > 0
  
  if (sum(keep.alleles) < 20) {next}
  
  c.mat.allele <- as.matrix(mat.allele[keep.alleles, ])
  c.mat.cov <- as.matrix(mat.cov[keep.alleles, ])
  
  #Get n of uniquely covered positions per sample pair
  mat.shared_cov[s, ] <- apply(c.mat.cov, 2, function(v) {length(unique(v))})
  
  #Sweep to get raw differences
  mat.diff[s, ] <- colSums(abs(sweep(c.mat.allele, 1, c.mat.allele[, s])) * (c.mat.cov > 0))
}

#Get sample combinations
combn.samples <- t(combn(data.sample$sample_alias, 2))

#Pull together results structure
tibble(
  species = species,
  sample_1 = combn.samples[, 1],
  sample_2 = combn.samples[, 2],
  shared_cov = as.dist(mat.shared_cov),
  diff_abs = as.dist(mat.diff),
  dist = as.dist(mat.diff) / as.dist(mat.shared_cov)
) %>%
  drop_na(shared_cov) %>%
  filter(dist < Inf & dist > -Inf) -> dat.allele_dist

#Store
save(dat.allele_dist, file = paste0(folder.base, "results/allele.dist/allele_distances.", species, ".RData"))

#Report & quit
writeLines(paste(date(), species, "=> DONE."))
q("no")
################################################################################
################################################################################







###
