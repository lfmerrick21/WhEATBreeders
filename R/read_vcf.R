# This function will read a VCF file into a vcf genotype table
read_vcf <- function(infile, skip = 8)
{
  vcf_table <- read.table(infile, header = T, skip = skip, stringsAsFactors = F, comment.char = "", na = ".")
  names(vcf_table)[1] <- "CHROM"

  return(vcf_table)
}
