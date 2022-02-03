# Perform mean imputation for a given genotype table
# genotypes
impute_mean <- function(genotypes)
{
  for (i in 1:ncol(genotypes))
  {
    na_indices <- is.na(genotypes[,i])
    genotypes[na_indices,i] <- mean(genotypes[,i], na.rm = T)
  }
  return(genotypes)
}
