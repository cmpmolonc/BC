# Stratified batch correction with PAM50 labels
# Samples have been previously assigned a PAM50 subtype. 
#
# @param exprs     Expression matrix (genes x samples), log2 normalised
# @param batch     Factor vector of batch (e.g., discovery vs validation)
# @param pam50     Factor vector of PAM50 labels, length = ncol(exprs)
# @return          Batch-corrected expression matrix
#
#
# Example usage: corrected <- stratified_batchcorrect(exprs, batch, pam50)

stratified_batchcorrect <- function(exprs, batch, pam50) {
  if (!requireNamespace("sva", quietly = TRUE)) {
    stop("Package 'sva' is required. Please install with install.packages('sva')")
  }
  
  # Ensure identical sample order
  if (!all(colnames(exprs) == names(batch))) {
    warning("Aligning sample order between exprs and batch/pam50")
    exprs <- exprs[, names(batch)]
  }
  
  # Create empty matrix to fill with correction values
  corrected <- matrix(NA, nrow = nrow(exprs), ncol = ncol(exprs))
  rownames(corrected) <- rownames(exprs)
  colnames(corrected) <- colnames(exprs)
  
  # Iterate over PAM50 subtype factor levels
  for (sub in levels(pam50)) {
    cat("Correcting subtype:", sub, "\n")
    
    idx <- which(pam50 == sub)
    if (length(unique(batch[idx])) < 2) {
      warning(paste("Subtype", sub, "only has one batch. Skipping correction."))
      corrected[, idx] <- exprs[, idx]
    } else {
      corrected[, idx] <- sva::ComBat(
        dat = exprs[, idx, drop = FALSE],
        batch = batch[idx],
        mod = NULL, # no covariates here
        par.prior = TRUE, #parametric correction
        prior.plots = FALSE # no plots
      )
    }
  }
  
  return(corrected)
}
