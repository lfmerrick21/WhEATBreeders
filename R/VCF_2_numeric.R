#' Convert VCF genotypes to GAPIT numeric format
#'
#' @param vcf_genotypes A dataframe of vcf formatted genotypes. The first 9 columns are marker metadata followed by genotypes. Columns = metadata and taxa, rows = markers
#'
#' @return A list containing a genotype dataframe and a marker map dataframe
VCF_2_numeric <- function(vcf_genotypes)
{
  meta_cols <- 9 # There are 9 marker data columns before the actual genotypes
  n_markers <- nrow(vcf_genotypes)
  n_samples <- ncol(vcf_genotypes)-meta_cols

  # Create map table
  marker_map <- data.frame(matrix(ncol = 3, nrow = nrow(vcf_genotypes)))
  names(marker_map) <- c("rs", "chr", "pos")
  # Extract map data
  marker_map$rs <- vcf_genotypes$ID
  marker_map$chr <- vcf_genotypes$CHROM
  marker_map$pos <- vcf_genotypes$POS

  # Extract vcf call data
  vcf_calls <- vcf_genotypes[,(meta_cols+1):ncol(vcf_genotypes)]
  row.names(vcf_calls) <- vcf_genotypes$ID

  t_genotypes <- matrix(ncol = ncol(vcf_calls), nrow = nrow(vcf_calls))
  # Convert unphased genotypes
  t_genotypes[vcf_calls == "0/0"] <- 0
  t_genotypes[vcf_calls == "0/1"] <- 1
  t_genotypes[vcf_calls == "1/0"] <- 1
  t_genotypes[vcf_calls == "1/1"] <- 2

  # Convert phased genotypes
  t_genotypes[vcf_calls == "0|0"] <- 0
  t_genotypes[vcf_calls == "0|1"] <- 1
  t_genotypes[vcf_calls == "1|0"] <- 1
  t_genotypes[vcf_calls == "1|1"] <- 2

  # Transpose the t_genotypes to match GAPIT specs
  num_genotypes <- data.frame(t(t_genotypes), stringsAsFactors = F)
  names(num_genotypes) <- row.names(vcf_calls)
  row.names(num_genotypes) <- names(vcf_calls)

  # Return a list of the genotypes and marker map
  convert_list <- list(num_genotypes, marker_map)
  names(convert_list) <- c('genotypes', 'marker_map')
  return(convert_list)
}
