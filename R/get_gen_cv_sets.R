# This function uses a cvFolds assignment to set matrix values to missing and retains the original values
# folds = a cvFolds object that split all genotype indices into different groups
# genotypes = a numeric genotype matrix (rows = taxa, columns = markers)
# fold_ct = the number of folds to use for genotype subsetting (if NULL, use all folds)
get_gen_cv_sets <- function(folds, genotypes, fold_ct = NULL)
{
  cv_genotype_list <- list()

  if (!is.null(fold_ct))
  {
    iterator <- 1:fold_ct
  } else
  {
    iterator <- unique(folds$which)
  }

  for (i in iterator)
  {
    fold_genotypes <- unlist(genotypes) # Decompose the matrix to simplify cross-validation index reference
    fold_genotypes[folds$subsets[folds$which == i]] <- NA # Set the values to missing
    fold_genotypes <- data.frame(matrix(fold_genotypes, ncol = ncol(genotypes), nrow = nrow(genotypes))) # Recompose the genotype matrix
    names(fold_genotypes) <- names(genotypes)
    row.names(fold_genotypes) <- row.names(genotypes)

    fold_list <- list(fold_genotypes)
    names(fold_list) <- c('fold_genotypes')

    cv_genotype_list[[i]] <- fold_list
    names(cv_genotype_list)[i] <- paste('fold_', i)
  }

  return(cv_genotype_list)
}
