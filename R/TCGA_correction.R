#   exprs:   gene x sample matrix (Variance-stabilised TCGA-RNASeq expression), rownames=genes, colnames=sample IDs
#   batch:   factor of batch (e.g., "1","2") corresponding to date of flow cell processing
#   pam50:   PAM50 sample assignment  (e.g., "LumA","LumB","Basal","Her2","Normal") 


exprs<-assay(TCGAVSD) # expression values of variance-stabilised TCGA-BRCA
batch<-readRDS("TCGA_batch_factor.rds") 
pam50<-readRDS("TCGA_PAM50_factor.rds")


stopifnot(all(colnames(exprs) %in% names(batch)))
exprs <- exprs[, names(batch)]  # ensure alignment of columns to rows exactly


####4 Combat Variations
## 1) ComBat: Parametric (no covariates) 
combat_parametric <- ComBat(
  dat         = exprs,
  batch       = batch,
  mod         = NULL,
  par.prior   = TRUE,    # parametric
  prior.plots = FALSE
)

## 2) ComBat: Non-parametric (no covariates)
combat_nonparam <- ComBat(
  dat         = exprs,
  batch       = batch,
  mod         = NULL,
  par.prior   = FALSE,   # non-parametric
  prior.plots = FALSE
)

## 3) ComBat: With PAM50 covariate
mod_pam50 <- model.matrix(~ pam50)  # intercept + PAM50 labels
combat_cov_pam50 <- ComBat(
  dat         = exprs,
  batch       = batch,
  mod         = mod_pam50,
  par.prior   = TRUE,
  prior.plots = FALSE
)

## 4) PAM50-stratified ComBat (per-subtype correction) 
## Use function defined in stratified_batchcorrect.R
combat_stratified <- stratified_batchcorrect(exprs, batch = batch, pam50 = pam50)

## Collect corrected outputs
combat_results <- list(
  parametric_noCov   = combat_parametric,
  nonparam_noCov     = combat_nonparam,
  covariate_PAM50    = combat_cov_pam50,
  stratified_byPAM50 = combat_stratified
)

# Access individual results as below
# exprs_covariate_corrected <- combat_results$covariate_PAM50
