# Calculates the minor allele frequency for every marker in a genotype matrix (coded as c(-1,0,1))
calc_maf_apply <- function(gt_mat, encoding = c(-1, 0, 1))
{
  col_func1 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == -1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)
    allele2_ct <- (sum(gt_col == 1, na.rm = T) * 2) + sum(gt_col == 0, na.rm = T)

    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }

  col_func2 <- function(gt_col)
  {
    allele1_ct <- (sum(gt_col == 0, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)
    allele2_ct <- (sum(gt_col == 2, na.rm = T) * 2) + sum(gt_col == 1, na.rm = T)

    maf <- min(c(allele1_ct, allele2_ct)) / (sum(!is.na(gt_col))*2)
  }

  if (all(encoding == c(-1, 0, 1)))
  {
    maf_vect <- apply(gt_mat, 2, col_func1)
  } else if (all(encoding == c(0, 1, 2)))
  {
    maf_vect <- apply(gt_mat, 2, col_func2)
  } else{
    print('Encoding not recognized, returning NULL')
    maf_vect <- NULL
  }

  return(maf_vect)
}
