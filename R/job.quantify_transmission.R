#/usr/bin/Rscript
################################################################################
#Parse and consolidate SNV data per pangenome
#
#=> process data for one species at a time
#=> load sample and profile data
#=> load gene content data
#=> load SNV data
#=> quantify engraftment based on gene content and SNVs
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
options(dplyr.summarise.inform = FALSE)
################################################################################
################################################################################


################################################################################
################################################################################
#Get current species ID
args <- commandArgs(TRUE)
species <- args[2]
mode <- args[3]

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
file.data_SNV.filtered <- paste0(folder.gene_snv, species, ".snv_filtered.Rdata")

#Check if the expected files exist, and quit if they don't
if (!file.exists(file.data_cov.filtered)) {writeLines(paste(date(), species, "=> filtered gene coverage file does not exist")); q("no")}
if (!file.exists(file.data_SNV.filtered)) {writeLines(paste(date(), species, "=> SNV file does not exist. Skipping.")); q("no")}

#Load input files
load(file.data_sample)
load(file.data_mOTU)
load(file.data_genome)
load(file.genome_sets)
load(file.data_cov.filtered)
load(file.data_SNV.filtered)

#Coerce data.fmt into data.frame
data.fmt <- as.data.frame(data.fmt)

#Extrapolate missing mOTU profile for one sample
data.mOTU <- cbind(data.mOTU, NA); colnames(data.mOTU)[ncol(data.mOTU)] <- "SAMN06018976"
data.mOTU.rel <- cbind(data.mOTU.rel, NA); colnames(data.mOTU.rel)[ncol(data.mOTU.rel)] <- "SAMN06018976"

#Fix erroneous entry for donor sample in data.fmt
#=> FMT_Smilie had replicates that were merged upstream; data.fmt may still point to the wrong ID
data.fmt[data.fmt$sample.donor == "SAMEA104393714" & !is.na(data.fmt$sample.donor), "sample.donor"] <- "SAMEA104393713"

#Calculate the total coverage per position and sample
#=> by summing allele coverages per position
dat.obs <- as.data.table(dat.snv %>% select(-base, -freq.global))[, lapply(.SD, sum, na.rm=T), by = c("gene", "pos", "pos.full")] %>% as_tibble()
gc(verbose = F)

#Get global prevalence of >=0 coverage per each position
dat.obs$prev.samples <- rowSums(dat.obs[, 4:ncol(dat.obs)] > 0)
################################################################################
################################################################################


################################################################################
################################################################################
#Quantify engraftment, rejection and co-existence based on SNVs
#=> get donor profile, pre-FMT recipient profile and post-FMT recipient profile
#=> define "determinant" SNVs (sensu Simone) and states (including confident absence of SNVs/genes)
#=> adjust SNV weight by "coreness" of genes in species?
#=> summarise SNV haplotypes per gene per sample
################################################################################
################################################################################
#Preallocate results collector
res.SNV <- data.frame()

#Loop through FMTs
#=> placebo FMTs will be treated differently, below
for (fmt in data.fmt$fmt.id) {
  #Skip if this is a placebo/autologous FMT
  if (data.fmt[fmt, "fmt.type"] == "placebo") {next()}
  
  if (mode == "background" & data.fmt[fmt, "fmt.type"] == "FMT") {next()}
  if (mode == "observation" & data.fmt[fmt, "fmt.type"] != "FMT") {next()}
  
  #Get current donor sample
  s.donor <- data.fmt[fmt, "sample.donor"]
  #Skip if donor sample was not available
  if (is.na(s.donor)) {next()}
  
  #Get curret pre-FMT sample
  s.pre <- data.fmt[fmt, "sample.pre_fmt"]
  #Skip if pre-FMT sample was not available
  if (is.na(s.pre)) {next()}
  
  #Subset sample data to current FMT
  if (data.fmt[fmt, "fmt.type"] == "FMT") {
    d.sample <- data.sample %>% filter(fmt.id == fmt | sample_alias %in% c(s.donor, s.pre))
  } else {
    d.sample <- data.sample %>% filter(sample_alias %in% samples.fmt[[fmt]])
  }
  
  #Double check for missing SNV data
  if (! all(d.sample$sample_alias %in% colnames(dat.obs))){
    samples.keep <- d.sample$sample_alias[d.sample$sample_alias %in% colnames(dat.obs)]
    samples.drop <- d.sample$sample_alias[! d.sample$sample_alias %in% colnames(dat.obs)]
    
    d.sample <- d.sample %>% filter(sample_alias %in% samples.keep)
    
    #Report
    for (dr.smpl in samples.drop) {writeLines(paste(date(), "=> drop sample", dr.smpl, "from FMT", fmt, "in", species))}
    
    #Skip current iteration if pre-FMT or donor sample were concerned
    if (s.pre %in% samples.drop | s.donor %in% samples.drop) {next()}
  }
  
  #Get abundance of current species across current samples
  if (species %in% rownames(data.mOTU.rel)) {
    d.sample$abundance <- data.mOTU.rel[species, d.sample$sample_alias]
  } else {
    d.sample$abundance <- colSums(dat.pres.gene[, d.sample$sample_alias])
    d.sample$abundance[d.sample$abundance != 0] <- NA
  }
  
  #Prepare and pre-process donor and recipient SNV profiles
  d.inc.gene <- dat.inc.gene[unique(dat.obs$gene), d.sample$sample_alias]
  d.pres.gene <- dat.pres.gene[unique(dat.obs$gene), d.sample$sample_alias]
  tmp.obs <- dat.obs[, c("gene", d.sample$sample_alias)] %>% group_by(gene) %>% summarise(across(.fns = sum))
  
  keep.genes <- rownames(d.inc.gene)[rowSums(d.inc.gene) > 0 & rowSums(tmp.obs[, 2:ncol(tmp.obs)]) > 1]
  
  #Skip if <3 genes have useful/sufficient coverage
  if (length(keep.genes) < 3) {next()}
  
  #Apply gene presence filter temporarily
  tmp.obs <- dat.obs[, c("gene", "pos", "pos.full", d.sample$sample_alias)] %>% filter(gene %in% keep.genes)
  d.inc.gene <- d.inc.gene[keep.genes, ]
  d.pres.gene <- d.pres.gene[keep.genes, ]
  
  #Remove all positions with <2x coverage across samples
  keep.pos <- tmp.obs$pos.full[rowSums(tmp.obs[, 4:ncol(tmp.obs)]) > 1]
  if (length(keep.pos) < 50) {next()}
  
  #Filter SNVs by gene presence and SNV observation
  #=> remove SNVs in genes not observed in any sample for the current FMT
  #=> remove positions with no coverage in any of the present samples
  d.snv <- dat.snv[, c("gene", "pos", "base", "pos.full", d.sample$sample_alias)] %>% filter(gene %in% keep.genes) %>% filter(pos.full %in% keep.pos)
  d.obs <- tmp.obs %>% filter(pos.full %in% keep.pos)
  d.inc.gene <- d.inc.gene[rownames(d.inc.gene) %in% unique(d.obs$gene), ]
  d.pres.gene <- d.pres.gene[rownames(d.pres.gene) %in% unique(d.obs$gene), ]
  
  #Iterate through samples and generate allele tables
  l.snv <- list()
  for (s in d.sample$sample_alias) {
    d.snv %>%
      select(gene, pos, base, pos.full, !!s) %>%
      rename(snv=s) %>%
      spread(base, snv, fill=0) %>%
      left_join(as_tibble(d.pres.gene, rownames = "gene") %>% select(gene, !!s) %>% rename(pres=s), by="gene") %>%
      mutate(cov=A+C+G+T) %>%
      #mutate(X=ifelse(pres == 0, 1, ifelse(pres == -1, NA, 0))) %>%
      mutate(X=ifelse(pres == 0 & cov < 1, 1, 0)) %>%
      mutate(
        A.tmp=ifelse(A > 0, "A", ""),
        C.tmp=ifelse(C > 0, "C", ""),
        G.tmp=ifelse(G > 0, "G", ""),
        T.tmp=ifelse(T > 0, "T", ""),
        X.tmp=ifelse(X == 1, "X", "")
      ) %>%
      mutate(SNV=str_c("", A.tmp, C.tmp, G.tmp, T.tmp, X.tmp)) %>%
      select(gene, pos, pos.full, cov, A, C, G, T, X, SNV) -> l.snv[[s]]
    
    #Override spurious observations in confidently absent genes
    l.snv[[s]]$SNV[is.na(l.snv[[s]]$SNV)] <- ""
    #l.snv[[s]]$SNV[grepl("X", l.snv[[s]]$SNV)] <- "X"
  }
  
  #Iterate through post-FMT samples and quantify engraftment, rejection and co-existence signals
  for (s in d.sample$sample_alias[d.sample$timepoint.fmt > 0]) {
    #Get current profiles
    p.donor <- l.snv[[s.donor]]
    p.pre <- l.snv[[s.pre]]
    p.post <- l.snv[[s]]
    
    #Identify "determinant" SNVs
    #=> first, exclude "trivial" positions, i.e. those where all three samples agree fully
    pos.trivial <- p.donor$SNV == p.pre$SNV & p.donor$SNV == p.post$SNV
    
    #=> get candidate determinant positions (coverage >=2 or confidently absent in all three samples)
    pos.candidate <- (p.donor$cov > 1 | p.donor$SNV == "X") & (p.pre$cov > 1 | p.pre$SNV == "X") & (p.post$cov > 1 | p.post$SNV == "X") & !pos.trivial
    if (sum(pos.candidate, na.rm=T) < 50) {next()}
    
    #Subset data
    p.donor <- p.donor[pos.candidate, ] %>% select(-SNV) %>% gather(snv, count, A, C, G, T) %>% mutate(count=ifelse(X > 0, 0, ifelse(count == 0, -1, count))) %>% select(-X) %>% arrange(gene, pos)
    p.pre <- p.pre[pos.candidate, ]  %>% select(-SNV) %>% gather(snv, count, A, C, G, T) %>% mutate(count=ifelse(X > 0, 0, ifelse(count == 0, -1, count))) %>% select(-X) %>% arrange(gene, pos)
    p.post <- p.post[pos.candidate, ]  %>% select(-SNV) %>% gather(snv, count, A, C, G, T) %>% mutate(count=ifelse(X > 0, 0, ifelse(count == 0, -1, count))) %>% select(-X) %>% arrange(gene, pos)
    
    #Re-filter to remove SNVs which are trivially absent in all three samples
    snv.candidate <- !(p.donor$count <= 0 & p.pre$count <= 0 & p.post$count <= 0)
    if (sum(snv.candidate, na.rm=T) < 50) {next()}
    p.donor <- p.donor[snv.candidate, ]
    p.pre <- p.pre[snv.candidate, ]
    p.post <- p.post[snv.candidate, ]
    
    #Identify determinant features
    #=> SNV or confident absence unique to donor, pre-FMT or post-FMT recipient samples
    #=> determinant donor: unique relative to pre-FMT
    #=> determinant pre-FMT: unique relative to donor
    #=> determinant post-FMT: unique relative to donor & pre-FMT
    det.donor <- p.donor$cov > 0 & p.donor$count > 0 & ((p.pre$cov > 0 & p.pre$count <= 0) | (p.pre$cov == 0 & p.pre$count == 0))
    det.pre <- p.pre$cov > 0 & p.pre$count > 0 & ((p.donor$cov > 0 & p.donor$count <= 0) | (p.donor$cov == 0 & p.donor$count == 0))
    det.post <- p.post$cov > 0 & p.post$count > 0 & ((p.donor$cov > 0 & p.donor$count <= 0) | (p.donor$cov == 0 & p.donor$count == 0)) & ((p.pre$cov > 0 & p.pre$count <= 0) | (p.pre$cov == 0 & p.pre$count == 0))
    
    if (sum(det.donor | det.pre | det.post) < 50) {next()}
    
    #Count events / positions supporting engraftment, rejection, co-existence and environmental influx patterns
    cbind(p.post, data.frame(det.donor=det.donor, count.donor=p.donor$count, det.pre=det.pre, count.pre=p.pre$count, det.post=det.post)) %>%
      filter(det.donor | det.pre | det.post) %>%
      group_by(gene, pos, pos.full, snv) %>%
      summarise(
        #Engraftment events: donor determinant position and SNV/feature observed in post-FMT sample
        #=> pattern [D - R/X - D] or [X - R - X]
        engraftment = ifelse((det.donor | det.pre) & ((count.donor > 0 & count.pre <= 0 & count > 0) | (count.donor == 0 & count.pre > 0 & count == 0)), 1, 0),
        #Rejection events: pre-FMT determinant positions and SNV/feature observed in post-FMT sample
        #=> pattern [D/X - R - R] or [D - X - X]
        rejection = ifelse((det.donor | det.pre) & ((count.pre > 0 & count.donor <= 0 & count > 0) | (count.pre == 0 & count.donor > 0 & count == 0)), 1, 0),
        #Novel signature events: post-FMT determinant positions
        influx = ifelse(det.post, 1, 0)
      ) %>%
      group_by(gene, pos, pos.full) %>%
      summarise(
        engraftment = sum(engraftment),
        rejection = sum(rejection),
        influx = sum(influx)
      ) %>%
      ungroup() %>%
      mutate(
        donor = ifelse(engraftment > 0 & rejection == 0, engraftment / (engraftment + influx), 0),
        recipient = ifelse(engraftment == 0 & rejection > 0, rejection / (rejection + influx), 0),
        coexistence = ifelse(engraftment > 0 & rejection > 0, (engraftment + rejection) / (engraftment + rejection + influx), 0),
        influx = ifelse(influx > 0, influx / (engraftment + rejection + influx), 0)
      ) %>%
      select(-engraftment, -rejection) %>%
      summarise(
        donor = sum(donor),
        recipient = sum(recipient),
        coexistence = sum(coexistence),
        influx = sum(influx)
      ) -> c.res
    
    #Pass results
    res.SNV <- rbind(
      res.SNV,
      data.frame(
        species = species,
        fmt.id = fmt,
        fmt.type = data.fmt[fmt, "fmt.type"],
        timepoint.fmt = d.sample[d.sample$sample_alias == s, "timepoint.fmt"],
        sample.donor = s.donor,
        sample.pre = s.pre,
        sample.post = s,
        abd.donor = d.sample[d.sample$sample_alias == s.donor, "abundance"],
        abd.pre = d.sample[d.sample$sample_alias == s.pre, "abundance"],
        abd.post = d.sample[d.sample$sample_alias == s, "abundance"],
        n.pos.trivial = sum(pos.trivial),
        n.pos.candidate = sum(pos.candidate),
        n.det.donor = sum(det.donor),
        n.det.pre = sum(det.pre),
        n.det.post = sum(det.post),
        engraftment = c.res$donor,
        rejection = c.res$recipient,
        coexistence = c.res$coexistence,
        influx = c.res$influx
      )
    )
  }
}
################################################################################
################################################################################
#Loop through "placebo" (autologous) FMTs
#=> these are treated differently
for (fmt in data.fmt[data.fmt$fmt.type == "placebo", "fmt.id"]) {
  if (mode == "background") {next()}
  
  #Get curret pre-FMT sample
  s.pre <- data.fmt[fmt, "sample.pre_fmt"]
  #Skip if pre-FMT sample was not available
  if (is.na(s.pre)) {next()}
  
  #Subset sample data to current FMT
  d.sample <- data.sample %>% filter(fmt.id == fmt | sample_alias %in% c(s.pre))
  
  #Double check for missing SNV data
  if (! all(d.sample$sample_alias %in% colnames(dat.obs))){
    samples.keep <- d.sample$sample_alias[d.sample$sample_alias %in% colnames(dat.obs)]
    samples.drop <- d.sample$sample_alias[! d.sample$sample_alias %in% colnames(dat.obs)]
    
    if (length(samples.keep) <= 1) {next()}
    
    d.sample <- d.sample %>% filter(sample_alias %in% samples.keep)
    
    #Report
    for (dr.smpl in samples.drop) {writeLines(paste(date(), "=> drop sample", dr.smpl, "from FMT", fmt, "in", species))}
    
    #Skip current iteration if pre-FMT or donor sample were concerned
    if (s.pre %in% samples.drop) {next()}
  }
  
  #Get abundance of current species across current samples
  if (species %in% rownames(data.mOTU.rel)) {
    d.sample$abundance <- data.mOTU.rel[species, d.sample$sample_alias]
  } else {
    d.sample$abundance <- colSums(dat.pres.gene[, d.sample$sample_alias])
    d.sample$abundance[d.sample$abundance != 0] <- NA
  }
  
  #Prepare and pre-process SNV profiles
  d.inc.gene <- dat.inc.gene[unique(dat.obs$gene), d.sample$sample_alias]
  d.pres.gene <- dat.pres.gene[unique(dat.obs$gene), d.sample$sample_alias]
  tmp.obs <- dat.obs[, c("gene", d.sample$sample_alias)] %>% group_by(gene) %>% summarise(across(.fns = sum))
  
  keep.genes <- rownames(d.inc.gene)[rowSums(d.inc.gene) > 0 & rowSums(tmp.obs[, 2:ncol(tmp.obs)]) > 1]
  
  #Skip if <3 genes have useful/sufficient coverage
  if (length(keep.genes) < 3) {next()}
  
  #Apply gene presence filter temporarily
  tmp.obs <- dat.obs[, c("gene", "pos", "pos.full", d.sample$sample_alias)] %>% filter(gene %in% keep.genes)
  d.inc.gene <- d.inc.gene[keep.genes, ]
  d.pres.gene <- d.pres.gene[keep.genes, ]
  
  #Remove all positions with <2x coverage across samples
  keep.pos <- tmp.obs$pos.full[rowSums(tmp.obs[, 4:ncol(tmp.obs)]) > 1]
  if (length(keep.pos) < 50) {next()}
  
  #Filter SNVs by gene presence and SNV observation
  #=> remove SNVs in genes not observed in any sample for the current FMT
  #=> remove positions with no coverage in any of the present samples
  d.snv <- dat.snv[, c("gene", "pos", "base", "pos.full", d.sample$sample_alias)] %>% filter(gene %in% keep.genes) %>% filter(pos.full %in% keep.pos)
  d.obs <- tmp.obs %>% filter(pos.full %in% keep.pos)
  d.inc.gene <- d.inc.gene[rownames(d.inc.gene) %in% unique(d.obs$gene), ]
  d.pres.gene <- d.pres.gene[rownames(d.pres.gene) %in% unique(d.obs$gene), ]
  
  #Iterate through samples and generate allele tables
  l.snv <- list()
  for (s in d.sample$sample_alias) {
    d.snv %>%
      select(gene, pos, base, pos.full, !!s) %>%
      rename(snv=s) %>%
      spread(base, snv, fill=0) %>%
      left_join(as_tibble(d.pres.gene, rownames = "gene") %>% select(gene, !!s) %>% rename(pres=s), by="gene") %>%
      mutate(cov=A+C+G+T) %>%
      #mutate(X=ifelse(pres == 0, 1, ifelse(pres == -1, NA, 0))) %>%
      mutate(X=ifelse(pres == 0 & cov < 1, 1, 0)) %>%
      mutate(
        A.tmp=ifelse(A > 0, "A", ""),
        C.tmp=ifelse(C > 0, "C", ""),
        G.tmp=ifelse(G > 0, "G", ""),
        T.tmp=ifelse(T > 0, "T", ""),
        X.tmp=ifelse(X == 1, "X", "")
      ) %>%
      mutate(SNV=str_c(A.tmp, C.tmp, G.tmp, T.tmp, X.tmp)) %>%
      select(gene, pos, pos.full, cov, A, C, G, T, X, SNV) -> l.snv[[s]]
    
    #Override spurious observations in confidently absent genes
    #l.snv[[s]]$SNV[grepl("X", l.snv[[s]]$SNV)] <- "X"
  }
  
  #Iterate through post-FMT samples and quantify engraftment, rejection and co-existence signals
  for (s in d.sample$sample_alias[d.sample$timepoint.fmt > 0]) {
    #Get current profiles
    p.pre <- l.snv[[s.pre]]
    p.post <- l.snv[[s]]
    
    #Identify "determinant" SNVs
    #=> first, exclude "trivial" positions, i.e. those where all samples agree fully
    pos.trivial <- p.pre$SNV == p.post$SNV
    
    #=> get candidate determinant positions (coverage >=2 or confidently absent in all three samples)
    pos.candidate <- (p.pre$cov > 1 | p.pre$SNV == "X") & (p.post$cov > 1 | p.post$SNV == "X") & !pos.trivial
    if (sum(pos.candidate) < 50) {next()}
    
    #Subset data
    p.pre <- p.pre[pos.candidate, ]  %>% select(-SNV) %>% gather(snv, count, A, C, G, T) %>% mutate(count=ifelse(X > 0, 0, ifelse(count == 0, -1, count))) %>% select(-X) %>% arrange(gene, pos)
    p.post <- p.post[pos.candidate, ]  %>% select(-SNV) %>% gather(snv, count, A, C, G, T) %>% mutate(count=ifelse(X > 0, 0, ifelse(count == 0, -1, count))) %>% select(-X) %>% arrange(gene, pos)
    
    #Re-filter to remove SNVs which are trivially absent in all three samples
    snv.candidate <- !(p.pre$count <= 0 & p.post$count <= 0)
    if (sum(snv.candidate) < 0) {next()}
    p.pre <- p.pre[snv.candidate, ]
    p.post <- p.post[snv.candidate, ]
    
    #Identify determinant features
    det.pre <- p.pre$cov > 0 & p.pre$count > 0
    det.post <- p.post$cov > 0 & p.post$count > 0 & ((p.pre$cov > 0 & p.pre$count <= 0) | (p.pre$cov == 0 & p.pre$count == 0))
    
    if (sum(det.pre | det.post) < 50) {next()}
    
    #Count events / positions supporting engraftment, rejection, co-existence and environmental influx patterns
    cbind(p.post, data.frame(det.pre=det.pre, count.pre=p.pre$count, det.post=det.post)) %>%
      filter(det.pre | det.post) %>%
      group_by(gene, pos, pos.full, snv) %>%
      summarise(
        #Engraftment events: n/a for placebo group
        engraftment = 0,
        #Rejection events: pre-FMT determinant positions and SNV/feeature observed in post-FMT sample
        rejection = ifelse(det.pre & ((count.pre > 0 & cov > 0 & count > 0) | (count.pre == 0 & count == 0)), 1, 0),
        #Novel signature events: post-FMT determinant positions
        influx = ifelse(det.post, 1, 0)
      ) %>%
      group_by(gene, pos, pos.full) %>%
      summarise(
        engraftment = sum(engraftment),
        rejection = sum(rejection),
        influx = sum(influx)
      ) %>%
      ungroup() %>%
      mutate(
        donor = ifelse(engraftment > 0 & rejection == 0, engraftment / (engraftment + influx), 0),
        recipient = ifelse(engraftment == 0 & rejection > 0, rejection / (rejection + influx), 0),
        coexistence = ifelse(engraftment > 0 & rejection > 0, (engraftment + rejection) / (engraftment + rejection + influx), 0),
        influx = ifelse(influx > 0, influx / (engraftment + rejection + influx), 0)
      ) %>%
      select(-engraftment, -rejection) %>%
      summarise(
        donor = sum(donor),
        recipient = sum(recipient),
        coexistence = sum(coexistence),
        influx = sum(influx)
      ) -> c.res
    
    #Pass results
    res.SNV <- rbind(
      res.SNV,
      data.frame(
        species = species,
        fmt.id = fmt,
        fmt.type = data.fmt[fmt, "fmt.type"],
        timepoint.fmt = d.sample[d.sample$sample_alias == s, "timepoint.fmt"],
        sample.donor = NA,
        sample.pre = s.pre,
        sample.post = s,
        abd.donor = NA,
        abd.pre = d.sample[d.sample$sample_alias == s.pre, "abundance"],
        abd.post = d.sample[d.sample$sample_alias == s, "abundance"],
        n.pos.trivial = sum(pos.trivial),
        n.pos.candidate = sum(pos.candidate),
        n.det.donor = 0,
        n.det.pre = sum(det.pre),
        n.det.post = sum(det.post),
        engraftment = c.res$donor,
        rejection = c.res$recipient,
        coexistence = c.res$coexistence,
        influx = c.res$influx
      )
    )
  }
}
################################################################################
################################################################################


################################################################################
################################################################################
#Save & exit
if (mode == "background") {
  save(res.SNV, data.fmt.shuffled, samples.fmt, file = paste0(folder.base, "results/transmission.snv/", species, ".transmission_snv.shuffled.Rdata"))
} else {
  save(res.SNV, file = paste0(folder.base, "results/transmission.snv/", species, ".transmission_snv.Rdata"))
}
writeLines(paste(date(), "=> done with", mode, species))
################################################################################
################################################################################








#

