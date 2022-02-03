# This function will test the accuracy of mean imputation across a set of folds
test_mean_impute_cv <- function(genotypes, genotypes_na_folds)
{
  accuracy_vect <- vector(mode = 'numeric', length = length(genotypes_na_folds))

  for (i in 1:length(accuracy_vect))
  {
    print(paste("Testing fold: ", i, sep = ''))
    genotypes_na <- genotypes_na_folds[[i]][[1]] # Extract the genotypes with missing values
    genotypes_imp <- impute_mean(genotypes_na) # Perform imputation
    accuracy_vect[i] <- calc_impute_acc(genotypes, genotypes_na, genotypes_imp) # Calculate an accuracy
  }

  return(accuracy_vect)
}
