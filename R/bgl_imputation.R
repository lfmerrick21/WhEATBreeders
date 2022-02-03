#create a beagle impuation function to easily create results over different missing rates
bgl_imputation <- function(data,map,missing_rate)
{
  X.raw=data#create dataframe that is just the marker data
  X=X.raw #create a copy of the dataframe to work on
  #Set variables to be used to randomize missing values
  n=nrow(X)
  m=ncol(X)
  dp=m*n #dp gives the total # of data points - use for % calculation later
  uv=runif(dp)#Converts each uniform vector into a series of TRUE/FALSE where each TRUE is a missing variable
  mr=missing_rate # missing rate (how much percentage of data is NA)
  missing=uv<mr # missing = T/F
  #Format indicator as matrix
  index.m=matrix(missing,n,m)
  #Set missing values as NA
  X[index.m]=NA
  # Convert the genotypes to a vcf
  na_file <- "genotypes_bgl"
  vcf_na <- numeric_2_VCF(X, map)
  write_vcf(vcf_na, outfile = na_file)
  # Assign parameters
  genotype_file = "genotypes_bgl.vcf"
  outfile = "genotypes_bgl_imp"

  # Define a system command
  command1_prefix <- "java -Xmx1g -jar beagle.25Nov19.28d.jar"
  command_args <- paste(" gt=", genotype_file, " out=", outfile, sep = "")
  command1 <- paste(command1_prefix, command_args)
  test <- system(command1, intern = T)
  # Run BEAGLE using the system function, this will produce a gzip .vcf file
  Xvcf=read.vcfR("genotypes_bgl_imp.vcf.gz")
  Xbgl=data.frame(Xvcf@fix,Xvcf@gt)
  genotypes_imp <- VCF_2_numeric(Xbgl)[[1]]
  # Calculate accuracy between original and imputed values
  accuracy.r=cor(X.raw[index.m],genotypes_imp[index.m], use = 'complete.obs')
  index.match=X.raw==genotypes_imp
  index.mm=index.match&index.m
  accuracy.m=length(X[index.mm])/length(X[index.m])
  accbgl=c(accuracy.r, accuracy.m)
  return(accbgl)
}
