#A function that will filter a genotype matrix based on maf and missingness
#Calculates the proportion of missing data for every marker in a genotype matrix (mising data is NA)
calc_missrate <- function(gt_mat)
{
  col_func <- function(gt_col)
  {
    missrate <- sum(is.na(gt_col)) / length(gt_col)
    return(missrate)
  }

  missrate_vect <- apply(gt_mat, 2, col_func)

  return(missrate_vect)
}
