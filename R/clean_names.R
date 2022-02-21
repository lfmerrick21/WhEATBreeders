clean_names <- function(x){
  x <- as.character(x)
  x<-toupper(x)
  x <- gsub(pattern = " ", replacement = "-",x = x, fixed = TRUE)
  x <- gsub(pattern = "--", replacement = "-",x = x, fixed = TRUE)
  return(x)
}
