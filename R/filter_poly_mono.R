filter_poly_mono <- function(data, missing_rate)
{
  M <- data
  n <- ncol(M)
  a <- nrow(M)
  #| Convert 1s to NA
  vars <- matrix(NA, 1, n)
  for (checkM in 1:n){                               #| Find markers who are monomorphic
    disone <- var(M[,checkM], na.rm = TRUE)
    vars[1,checkM] <- disone

  }
  length(which(vars == 0))
  rm_mono <- which(vars == 0)
  M <- M[,-rm_mono]
  #remove missing markers
  n <- ncol(M)
  a <- nrow(M)
  TooManyBlanks <- NULL
  for (i in 1:n){
    m <- sum(M[,i]==1)
    if (m > a*missing_rate) {TooManyBlanks <- c(TooManyBlanks, i)}
  }
  M <- M[,-TooManyBlanks]                         #| Remove MARKERS with >20% Missing Data
  n <- ncol(M)
  a <- nrow(M)
  M[M==1] <- NA

  filt=M
  return(filt)
}
