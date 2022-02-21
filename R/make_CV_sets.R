# This is a function that will split data into a list of k-folds
make_CV_sets <- function(list_length, k = 5){
  rand_values <- rnorm(list_length)
  k_quantiles <- quantile(rand_values, 0:k/k)
  k_assign <- cut(rand_values, k_quantiles, include.lowest = T, labels = F)

  cv_list <- list()
  for (i in 1:k)
  {
    fold_assignment <- k_assign != i
    cv_list[[i]] <- fold_assignment
  }
  return(cv_list)
}
