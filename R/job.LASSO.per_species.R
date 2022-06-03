#/usr/bin/Rscript
################################################################################
#LASSO regressions per-species to predict FMT outcomes
#
#2022-03
#sebastian.schmidt@embl.de
################################################################################


################################################################################
################################################################################
# Load Packages
suppressPackageStartupMessages(library("tidyverse", warn.conflicts = F, quietly=T))
library("glmnet", warn.conflicts = F, quietly=T)
library("caret", warn.conflicts = F, quietly=T)
library("vip", warn.conflicts = F, quietly=T)
library("pROC", warn.conflicts = F, quietly = T)

#Make R behave marginally less moronic
options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
################################################################################
################################################################################


################################################################################
################################################################################
#Get current species ID
args <- commandArgs(TRUE)
sp <- args[2]

#Set folder names
PARAM <- list()

#Set relevant folder names
load("[...]/data/param.analysis_transmission.Rdata")
PARAM$folder.base <- "~/" #=> SET TO APPROPRIATE REPO FOLDER

#Load preprocessed data
load(paste0(PARAM$folder.base, "data/data.sample.RData"))
load(paste0(PARAM$folder.base, "data/data.phenotype.Rdata"))
load(paste0(PARAM$folder.base, "data/data.transmission_snv.Rdata"))
load(paste0(PARAM$folder.base, "data/data.transmission_snv.scored.Rdata"))
load(paste0(PARAM$folder.base, "data/data.transmission_snv.scored.consolidated.Rdata"))
load(paste0(PARAM$folder.base, "data/data.mOTU.RData"))
load(paste0(PARAM$folder.base, "data/data.coverage.COG.MAG.Rdata"))
load(paste0(PARAM$folder.base, "data/data.tree.Rdata"))
load(paste0(PARAM$folder.base, "data/data.beta_div.RData"))
load(paste0(PARAM$folder.base, "data/data.allele_distances.FMT.Rdata"))
load(paste0(PARAM$folder.base, "data/data.allele_diversity.Rdata"))
load(paste0(PARAM$folder.base, "data/data.gmm.gmgc.Rdata"))
load(paste0(PARAM$folder.base, "data/data.gmm.by_species.Rdata"))
load(paste0(PARAM$folder.base, "data/data.gmgc.KO.Rdata"))
load(paste0(PARAM$folder.base, "data/data.sample.RData"))
################################################################################
################################################################################



################################################################################
################################################################################
#Precompute KEGG KO profiles
data.ko.gmgc %>%
  dplyr::select(ko, sample_alias, count_norm) %>%
  dplyr::mutate(count_norm = count_norm + 10^-9) %>%
  tidyr::complete(ko, sample_alias, fill = list(count_norm = 10^-9)) %>%
  dplyr::mutate(
    count.clr = log2(count_norm / mean(count_norm, na.rm=T))
  ) %>%
  group_by(ko) %>%
  dplyr::mutate(
    count.z = (count_norm - mean(count_norm, na.rm = T)) / sd(count_norm, na.rm=T)
  ) %>%
  ungroup() %>%
  dplyr::select(sample_alias, ko, count.clr) %>%
  pivot_wider(names_from = "ko", values_from = "count.clr") -> data.ko.gmgc.norm.wide
################################################################################
################################################################################


################################################################################
################################################################################
#Function from
#https://rdrr.io/github/HuntsmanCancerInstitute/hciR/src/R/as_matrix.R
as_matrix <- function(x){
  if(!tibble::is_tibble(x) ) stop("x must be a tibble")
  y <- as.matrix.data.frame(x[,-1])
  rownames(y) <- x[[1]]
  y
}

#Prepare species abundance data
data.mOTU.clr <- log2((data.mOTU.rel + 10^-6) / mean(data.mOTU.rel + 10^-6))
data.mOTU.clr.subset_50 <- data.mOTU.clr[use.species.engraftment.gt_50[use.species.engraftment.gt_50 %in% rownames(data.mOTU.clr)], ] %>%
  t() %>%
  as_tibble(rownames = "sample")

#Prepare engraftment data
dat.transmission.subset_gt_50 <- dat.transmission.scored.engraftment %>%
  filter(!sample.post == "SAMEA7082346") %>%
  group_by(fmt.id) %>% top_n(1, timepoint.fmt) %>% ungroup() %>%
  filter(species %in% use.species.engraftment.gt_50) %>%
  dplyr::select(species, fmt.id, frac.don) %>%
  pivot_wider(values_from = frac.don, names_from = species, values_fill = 0)

#Define variables
variables.metabolism.pre <- c(
  "proteolysis.pre",
  "saccharolysis.pre",
  "lipolysis.pre",
  "mucin_degradation.pre",
  "propionate_production.pre",
  "butyrate_production.pre"
)

variables.metabolism.donor <- c(
  "proteolysis.donor",
  "saccharolysis.donor",
  "lipolysis.donor",
  "mucin_degradation.donor",
  "propionate_production.donor",
  "butyrate_production.donor"
)

variables.metabolism.post <- c(
  "proteolysis.post",
  "saccharolysis.post",
  "lipolysis.post",
  "mucin_degradation.post",
  "propionate_production.post",
  "butyrate_production.post"
)

#Set variable types
var_lasso <- list(
  "ex_ante.technical" = c("min_depth", "indication.rcdi", "indication.ibd", "indication.ibs", "indication.mets", "indication.mel", "indication.esbl", "fmt_sample_prep.fresh", "fmt_route.nasal", "bowel.prep.summary", "antibiotics.summary"),
  "ex_ante.community.tax" = c("richness.donor.scaled", "richness.pre.scaled", "fr.donor.scaled", "fr.pre.scaled", "dissimilarity.donor_pre"),
  "ex_ante.community.metabolic" = c("dist_gmm.donor_pre", variables.metabolism.pre, variables.metabolism.donor),
  "ex_ante.focal_species" = c("abd_ratio.donor_pre", "allele_dist.donor_pre", "allele_div.pre", "allele_div.donor", "phylo.complementarity", "ko.complementarity"),
  "post_hoc.community.tax" = c("richness.post.scaled", "richness.ratio.post_pre", "fr.post.scaled", "dissimilarity.post_pre"),
  "post_hoc.community.metabolic" = c("dist_gmm.post_pre", variables.metabolism.post)
)
################################################################################
################################################################################


################################################################################
################################################################################
#Define function to run LASSO regression and return eval results
run.lasso <- function(
  c.vars,
  c.response = "engraftment",
  set.train = c.train,
  set.test = c.test,
  family = "gaussian",
  var.name,
  fold = fold,
  species = sp,
  c.res = c.res
) {
  #Get data
  c.train.X <- set.train %>% dplyr::select(c.vars) %>% mutate_all(~replace_na(., min(., na.rm=T))) %>% as.matrix.data.frame()
  c.test.X <- set.test %>% dplyr::select(c.vars) %>% mutate_all(~replace_na(., min(., na.rm=T))) %>% as.matrix.data.frame()
  resp.train <- set.train %>% pull(c.response)
  resp.test <- set.test %>% pull(c.response)
  
  #For binary variables, check for appropriate spread in training set
  check.cols <- names(which(apply(c.train.X, 2, function(x) {all(x %in% c(0,1))})))
  if (length(check.cols) > 0) {
    tmp.spread <- colSums(c.train.X[, check.cols, drop = F] > 0)
    drop.cols <- names(which(tmp.spread < 0.2 * nrow(c.train.X) | tmp.spread > 0.8 * nrow(c.train.X)))
    
    c.vars <- c.vars[!c.vars %in% drop.cols]
    if (length(c.vars) < 2) {return(c.res)}
    
    c.train.X <- c.train.X[, !colnames(c.train.X) %in% drop.cols]
    c.test.X <- c.test.X[, !colnames(c.test.X) %in% drop.cols]
  }
  
  #Set lower and upper bounds for coefficient scan according to variable type
  lb <- -Inf
  ub <- Inf
  if (str_detect(var.name, "facilitation") | str_detect(var.name, "\\.pos")) {lb <- 0}
  if (str_detect(var.name, "exclusion") | str_detect(var.name, "\\.neg")) {ub <- 0}
  
  #Run cv'ed LASSO regression
  c.lasso <- cv.glmnet(c.train.X, resp.train, alpha = 1, family = family, upper.limits = ub, lower.limits = lb)
  
  #Get standardised coefficients at lambda.1se using the Agresti method
  sds <- apply(c.train.X, 2, sd)
  cs <- as.matrix(coef(c.lasso, s = "lambda.1se"))
  std_coefs <- cs[-1, 1] * sds
  
  #Prepare output
  c.res[["coef"]] <- bind_rows(c.res[["coef"]], tibble(
    species = species,
    resp.var = c.response,
    cv.fold = fold,
    family = family,
    var.type = var.name,
    term = c.vars,
    coef = std_coefs
  ))
  
  #Evaluate model fit and pass outputs
  c.predict.test <- as.numeric(predict(c.lasso, c.test.X, "lambda.1se"))
  c.res[["predict"]] <- bind_rows(c.res[["predict"]], tibble(
    species = species,
    resp.var = c.response,
    cv.fold = fold,
    family = family,
    var.type = var.name,
    response = resp.test,
    predicted = c.predict.test
  ))
  
  if (family == "binomial") {
    #Pass AUC
    c.res[["metric"]] <- bind_rows(c.res[["metric"]], tibble(
      species = species,
      resp.var = c.response,
      cv.fold = fold,
      family = family,
      var.type = var.name,
      metric = "auc",
      non.trivial = var(c.predict.test) != 0,
      value = metric_auc(resp.test, c.predict.test)
    ))
    
    #Pass ROC curve
    tmp.roc <- invisible(roc(response = resp.test, predictor = as.numeric(c.predict.test), ret = "coords"))
    c.res[["roc"]] <- bind_rows(c.res[["roc"]], tibble(
      species = species,
      resp.var = c.response,
      cv.fold = fold,
      family = family,
      var.type = var.name,
      sensitivity = tmp.roc$sensitivities,
      specificity = tmp.roc$specificities
    ))
    
    #Pass
    return(c.res)
  } else {
    #Pass metric
    c.res[["metric"]] <- bind_rows(c.res[["metric"]], tibble(
      species = species,
      resp.var = c.response,
      cv.fold = fold,
      family = family,
      var.type = var.name,
      metric = c("rho", "rsq", "mse", "rmse"),
      non.trivial = var(c.predict.test) != 0,
      value = c(
        cor(resp.test, c.predict.test),
        metric_rsquared(resp.test, c.predict.test),
        metric_mse(resp.test, c.predict.test),
        metric_rmse(resp.test, c.predict.test)
      )
    ))
    
    #Pass
    return(c.res)
  }
}
################################################################################
################################################################################


################################################################################
################################################################################
#Pre-allocate cross-validation folds
#=> 5-fold, 80:20
n.fold <- 5

#Minimum donor takeover threshold to consider a species
thresh.takeover <- 0.2

#Preallocate results collectors
collect.res.engraft.lasso.per_species <- list()
collect.res.engraft.lasso.metric <- collect.res.engraft.lasso.coef <- tibble()

data.phenotype.gmm %>% filter(mOTU == sp) %>% pull(label) -> c.label
writeLines(paste(date(), "=>", c.label))


##############################
#Get current data
c.dat <- dat.transmission.scored.outcome %>%
  filter(!sample.post == "SAMEA7082346") %>%
  group_by(fmt.id) %>% top_n(1, timepoint.fmt) %>% ungroup() %>%
  #=> subset data to current species
  filter(species == sp) %>%
  #=> select variables of interest
  dplyr::select(species, fmt.id, sample.donor, sample.pre, sample.post, inc.donor, inc.pre, inc.post, engraftment:oc.rec_turnover, unique(unlist(var_lasso))) %>%
  replace_na(list(allele_div.pre = 0, allele_div.donor = 0)) %>%
  #=> add species abundances in donor
  left_join(data.mOTU.clr.subset_50, by = c("sample.donor" = "sample")) %>%
  #=> add species abundances in recipient pre-FMT
  left_join(data.mOTU.clr.subset_50, by = c("sample.pre" = "sample"), suffix = c("", ".abd_pre")) %>%
  #=> add species abundances in recipient post-FMT
  left_join(data.mOTU.clr.subset_50, by = c("sample.post" = "sample"), suffix = c(".abd_donor", ".abd_post")) %>%
  #=> add engraftment per species (excluding the focal species)
  left_join(dat.transmission.subset_gt_50 %>% dplyr::select(-sp), by = "fmt.id") %>%
  #=> add full KO profiles (clr-normalised)
  #xleft_join(data.ko.gmgc.norm.wide, by = c("sample.pre" = "sample")) %>%
  #=> remove superfluous columns
  dplyr::select(-species, -fmt.id, -sample.donor, -sample.pre, -sample.post)

#Get total fraction of FMTs per outcome category for this species
c.dat %>%
  dplyr::summarise(
    oc.donor_takeover = sum(oc.donor_takeover) / sum(inc.donor),
    oc.donor_engraftment = sum(oc.donor_engraftment) / sum(inc.donor),
    oc.donor_rejection = sum(oc.donor_rejection) / sum(inc.donor),
    oc.rec_persistence = sum(oc.rec_persistence) / sum(inc.pre),
    oc.rec_resilience = sum(oc.rec_resilience) / sum(inc.pre),
    oc.rec_species_loss = sum(oc.rec_species_loss) / sum(inc.pre),
    oc.rec_turnover = sum(oc.rec_turnover) / sum(inc.post)
  ) -> c.frac
##############################

##############################
#Drop abundance-based variable for non-mOTU species (where abundance cannot be properly estimated)
if (grepl("ANI_AL", sp)) {species.variables <- var_lasso$ex_ante.focal_species[var_lasso$ex_ante.focal_species != "abd_ratio.donor_pre"]} else {species.variables <- var_lasso$ex_ante.focal_species}
#Define variable names for species abundances in donor, pre and engraftment submatrices
var.abd_donor <- c.dat %>% dplyr::select(ends_with("abd_donor")) %>% names()
var.abd_pre <- c.dat %>% dplyr::select(ends_with("abd_pre")) %>% names()
var.abd_post <- c.dat %>% dplyr::select(ends_with("abd_post")) %>% names()
#var.KO_pre <- c.dat %>% dplyr::select(starts_with("K", ignore.case = F)) %>% names()
var.engraftment <- c.dat %>%
  dplyr::select((starts_with("ref_") | starts_with("meta_") | starts_with("ANI")) & !ends_with("abd_donor") & !ends_with("abd_pre") & !ends_with("abd_post") &!contains(sp)) %>%
  names()

#Get list of relevant variables
list.var <- c(
  var_lasso[names(var_lasso) != "ex_ante.focal_species"],
  list(
    ex_ante.focal_species = species.variables,
    ex_ante.abd_donor = var.abd_donor,
    ex_ante.abd_donor.facilitation = var.abd_donor,
    ex_ante.abd_donor.exclusion = var.abd_donor,
    ex_ante.abd_pre = var.abd_pre,
    ex_ante.abd_pre.facilitation = var.abd_pre,
    ex_ante.abd_pre.exclusion = var.abd_pre,
    ex_ante.combined = c(
      var_lasso$ex_ante.technical,
      var_lasso$ex_ante.community.tax,
      var_lasso$ex_ante.community.metabolic,
      species.variables,
      var.abd_pre,
      var.abd_donor
    ),
    post_hoc.abd_post = var.abd_post,
    post_hoc.abd_post.pos = var.abd_post,
    post_hoc.abd_post.neg = var.abd_post,
    post_hoc.engraftment = var.engraftment,
    post_hoc.combined = c(
      var_lasso$post_hoc.community.tax,
      var_lasso$post_hoc.community.metabolic,
      var.abd_post,
      var.engraftment
    )
  )
)
##############################

##############################
#Preallocate local results collectors across cv folds
c.res <- list()
c.res$predict <- c.res$coef <- c.res$metric <- c.res$roc <- tibble()
#Iterate
for (k in 1:5) {
  #Iterate through variables (variable types)
  for (resp in c(
    "frac.don",
    "frac.rec",
    "frac.coexist",
    "oc.donor_takeover",
    "oc.donor_engraftment",
    "oc.donor_rejection",
    "oc.rec_persistence",
    "oc.rec_resilience",
    "oc.rec_species_loss",
    "oc.rec_turnover"
  )) {
    #Set model type based on response variable type
    if (str_detect(resp, "frac")) {
      family <- "gaussian"
    } else {
      family <- "binomial"
    }
    
    #Skip if outcome data is too unbalanced
    if (str_detect(resp, "oc\\.") && (c.frac[[resp]] > 0.8 | c.frac[[resp]] < 0.2)) {next}
    
    #Sub-select appropriate ref data for current variable
    #=> for example, engraftment-related variables should only be checked in FMTs where the species was present in the donor, etc.
    if (resp %in% c("frac.don", "oc.donor_takeover", "oc.donor_engraftment", "oc.donor_rejection")) {
      c.dat %>% filter(inc.donor) -> m.dat
    } else if (resp %in% c("frac.rec", "oc.rec_persistence", "oc.rec_resilience", "oc.rec_species_loss")) {
      c.dat %>% filter(inc.pre) -> m.dat
    } else {
      c.dat %>% filter(inc.post) -> m.dat
    }
    
    #Skip if too little data is available
    if(nrow(m.dat) < 40) {next}
    
    #Roll the dice, stratified by target response variable
    if (str_detect(resp, "frac")) {
      cv.folds <- createFolds(m.dat[[resp]], k = 5)
    } else {
      cv.folds <- createFolds(as.factor(m.dat[[resp]]), k = 5)
    }
    
    #Iterate through cross-validation folds
    for (fold in 1:n.fold) {
      ##############################
      #Model engraftment success
      #=> attempt to predict the fraction of donor strain in the recipient post FMT
      #=> using R^2 and RMSE as metrics (attempting to predict a continuous variable)
      #=> attempt to predict in which FMTs the species will successfully take over the community (binary outcome)
      #=> using AUC as metric
      #=> but only for species where the total fraction of successful engraftment is ≥15%
      ##############################
      #Get current training and test sets
      c.test <- m.dat[cv.folds[[fold]], ]
      c.train <- m.dat[unlist(cv.folds[1:n.fold != fold]), ]
      
      #Skip for trivial outcomes
      if (sum(c.test[[resp]]) == 0 | sum(c.test[[resp]]) == nrow(c.test)) {next}
      
      #Iterate through variable types and run LASSO
      for (var.name in names(list.var)) {
        try(c.res <- run.lasso(
          c.vars = list.var[[var.name]],
          c.response = resp,
          var.name = var.name,
          set.train = c.train,
          set.test = c.test,
          fold = (k-1) * n.fold + fold,
          family = family,
          species = sp,
          c.res = c.res
        ))
      }
    }
  }
}

#Store data
save(c.res, file = paste0(PARAM$folder.base, "results/lasso.per_species/", sp, ".lasso.Rdata"))
################################################################################
################################################################################


################################################################################








