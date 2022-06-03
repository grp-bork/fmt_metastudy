#/usr/bin/Rscript
################################################################################
#Parse and consolidate SNV data per pangenome
#
#=> process data for one species at a time
#=> load sample and profile data
#=> load gene content data
#=> load SNV data
#
#2019-07-19
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
file.data_mOTU <- paste0(folder.data, "data.mOTU.RData")
file.data_genome <- paste0(folder.data, "genome_data.combined.Rdata")
file.genome_sets <- paste0(folder.data, "data.genome_sets.Rdata")
file.data_cov.filtered <- paste0(folder.gene_cov, species, ".cov.filtered.RData")
file.data_SNV.raw <- paste0(folder.gene_snv, species, ".snv_raw.tsv.gz")
file.bam_list <- paste0(folder.base, "parameters/pangenome_SNV.", species, ".bam_list.txt")

#Check if the expected files exist, and quit if they don't
if (!file.exists(file.data_cov.filtered)) {writeLines(paste(date(), species, "=> filtered gene coverage file does not exist")); q("no")}
if (!file.exists(file.data_SNV.raw)) {writeLines(paste(date(), species, "=> raw SNV file does not exist. Skipping.")); q("no")}

#Load input files
load(file.data_sample)
load(file.data_mOTU)
load(file.data_genome)
load(file.genome_sets)
load(file.data_cov.filtered)

#Extrapolate missing mOTU profile for one sample
data.mOTU <- cbind(data.mOTU, NA); colnames(data.mOTU)[ncol(data.mOTU)] <- "SAMN06018976"
data.mOTU.rel <- cbind(data.mOTU.rel, NA); colnames(data.mOTU.rel)[ncol(data.mOTU.rel)] <- "SAMN06018976"

#Fix erroneous entry for donor sample in data.fmt
#=> FMT_Smilie had replicates that were merged upstream; data.fmt may still point to the wrong ID
data.fmt[data.fmt$sample.donor == "SAMEA104393714" & !is.na(data.fmt$sample.donor), "sample.donor"] <- "SAMEA104393713"
################################################################################
################################################################################


################################################################################
################################################################################
#Parse raw SNV data
snv.raw <- fread(file=file.data_SNV.raw, sep="\t")

#Read bam list
bam.list <- readLines(file.bam_list) %>%
  gsub(".+/", "", .) %>%
  gsub(paste0(species, "."), "", .) %>%
  gsub("\\.screened.+", "", .) %>%
  gsub("_16s.+", "", .) %>%
  gsub("_17s.+", "", .) %>%
  gsub("\\.sorted.+", "", .) %>%
  gsub("\\.unique", "", .)

#Which samples should be kept?
#=> erroneously, 
keep.bam <- bam.list %in% data.sample$sample_alias

#Process SNV labels
tmp.allele <- str_split_fixed(snv.raw[[1]], ":", 3)

#Filter to only include genes from the pre-filtered list
keep.alleles.by_gene <- tmp.allele[,1] %in% rownames(dat.gene)

if (sum(keep.alleles.by_gene) < 200) {writeLines(paste(date(), species, "=> too few alleles left after gene filter. Skipping.")); q("no")}

#Apply first round of filters
#=> samples to keep
#=> alleles in relevant genes
dat.snv <- snv.raw[keep.alleles.by_gene, which(keep.bam) +1, with=FALSE]
colnames(dat.snv) <- bam.list[keep.bam]
rm(snv.raw); invisible(gc(verbose = F))
dat.allele <- as.data.frame(tmp.allele[keep.alleles.by_gene, ])
dat.allele <- cbind(dat.allele, paste(dat.allele[, 1], dat.allele[, 2], sep=":"))
colnames(dat.allele) <- c("gene", "pos", "base", "pos.full")
rm(tmp.allele); invisible(gc(verbose = F))

#Remove all SNVs which are not supported by >= 2 reads in at least 2 samples
if (species %in% c("ref_mOTU_v2_0466", "ref_mOTU_v2_0898", "ref_mOTU_v2_0899")) {
  keep.alleles.by_reads <- rowSums(dat.snv > 1) > 3
} else {
  keep.alleles.by_reads <- rowSums(dat.snv > 1) > 1
}
tmp.allele.pos.table <- table(dat.allele[keep.alleles.by_reads, "pos.full"])
keep.pos.table <- names(tmp.allele.pos.table)[tmp.allele.pos.table > 1]
keep.pos.by_reads <- dat.allele[, "pos.full"] %in% keep.pos.table
if (sum(keep.pos.by_reads) < 100) {writeLines(paste(date(), species, "=> too few positions with >=2 reads in >=2 samples. Skipping.")); q("no")}
#Apply filter
dat.snv <- dat.snv[keep.pos.by_reads, ]
dat.allele <- dat.allele[keep.pos.by_reads, ]
rm(tmp.allele.pos.table); invisible(gc(verbose = F))

#Calculate global background frequency of each allele
dat.allele <- as.data.frame(dat.allele)
dat.allele$freq.global <- rowSums(dat.snv > 0)

#Re-merge data structures
dat.snv <- cbind(dat.allele, dat.snv)
dat.snv$pos <- as.numeric(dat.snv$pos)
invisible(gc(verbose = F))

#Save SNV data for later (re-)use
save(dat.snv, file = paste0(folder.gene_snv, species, ".snv_filtered.Rdata"))

#Report
writeLines(paste(date(), "=> done with", species))
################################################################################
################################################################################

