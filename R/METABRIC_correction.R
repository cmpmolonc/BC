#   exprs:   gene x sample matrix (log2-normalised Illumina HT-12 METABRIC expression), rownames=genes, colnames=sample IDs
#   pheno:   data.frame with rows = samples (rownames = sample IDs), columns include:
#            pheno$Batch  (e.g., "Discovery"/"Validation")
#            pheno$PAM50  (e.g., "LumA","LumB","Basal","Her2","Normal")

stopifnot(all(colnames(exprs) %in% rownames(pheno)))
exprs <- exprs[, rownames(pheno)]  # align columns to pheno rows exactly

# Factors
batch <- factor(pheno$Batch, levels = c("Discovery","Validation"))
names(batch) <- rownames(pheno)

pam50 <- factor(pheno$PAM50, levels = c("LumA","LumB","Basal","Her2","Normal"))
names(pam50) <- rownames(pheno)

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
## Use function defined in stratified_batccorrect.R
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
