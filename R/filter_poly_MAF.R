filter_poly_MAF <- function(data, missing_rate)
{
  M <- data
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
  #| Convert 1s to NA
  vars <- matrix(NA, 1, n)
  for (checkM in 1:n){                               #| Find markers who are monomorphic
    disone <- var(M[,checkM], na.rm = TRUE)
    vars[1,checkM] <- disone

  }
  length(which(vars == 0))
  rm_mono <- which(vars == 0)
  M <- M[,-rm_mono]
  K = str(rm_mono)
  filt=M
  minor_allele_freq <- apply(filt, 2, function(x)
  {
    allele_freq1 <- (sum(x == 0)*2 + sum(x == 1)) / (sum(!is.na(x)) * 2)
    allele_freq2 <- (sum(x == 2)*2 + sum(x == 1)) / (sum(!is.na(x)) * 2)
    return(min(allele_freq1, allele_freq2))
  })
  filt_MAF <- filt[,minor_allele_freq > 0.05]
  return(filt_MAF)
}
