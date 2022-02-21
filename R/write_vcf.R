# This function will write a VCF to a text file in the proper format
write_vcf <- function(vcf_genotypes, outfile)
{
  # Set constants
  VCF_version = 4.2 # (This might be incorporated as a feature in the future for different methods for different versions)

  sink(paste(outfile, ".vcf", sep = "")) # Redirect all R output to the specified file
  cat(paste("##fileformat=VCFv", VCF_version, "\n", sep = '')) # Write the first line representing theversion
  cat('#')
  sink() # Close redirection

  write.table(vcf_genotypes, file = paste(outfile, ".vcf", sep = ""), append = T, quote = F, sep = "\t", na = ".", row.names = F, col.names = T) # Append the vcf_table to the output file with proper seperators
}
