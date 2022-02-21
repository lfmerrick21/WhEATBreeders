boxcox_t=function(vector){
  library(forecast)
  # to find optimal lambda
  lambda = BoxCox.lambda(vector )
  # now to transform vector
  T_box = BoxCox(vector, lambda)
  return(T_box)
}
