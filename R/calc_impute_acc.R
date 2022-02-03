# Calculate imptuation accuracy
calc_impute_acc <- function(genotypes, genotypes_na, genotypes_imp)
{
  na_indices <- which(is.na(genotypes_na), arr.ind = T) # Get the indices of missing values

  orig_values <- genotypes[na_indices]
  imp_values <- genotypes_imp[na_indices]

  accuracy <- cor(orig_values, imp_values)

  return(accuracy)
}
