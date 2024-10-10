#' Compute nucleotide diversity for each sample
#'
#' This function computes the nucleotide diversity for each sample in a haploid dataset. It converts the provided genetic dataframe into a `DNAbin` object for each sample, calculates the nucleotide diversity, and returns the results.
#'
#' @param arp_df A dataframe where each row corresponds to a sequence. The dataframe must contain at least the following columns:
#'   - `title` : sample names or identifiers.
#'   - columns two and following : genetic data in haploid format.
#' @return A named list where each element contains the nucleotide diversity for one sample. The names of the list elements correspond to the unique sample identifiers in the `title` column of the input data frame.
#'
#' @details
#' The function performs the following steps:
#' 1. For each unique sample (identified by the `title` column), it subsets the dataframe to include only the rows corresponding to that sample.
#' 2. Converts the subsetted dataframe to a `DNAbin` object using the `as.DNAbin` function.
#' 3. Calculates nucleotide diversity for each `DNAbin` object using the `nuc.div` function.
#'
#' @examples
#' nuc_div_results <- haploid.nuc.div(arp_df)
#' print(nuc_div_results)
#'
#' @export
haploid.nuc.div <- function(arp_df)
{
  dnabin_list <- setNames(lapply(unique(arp_df$title), function(pop)
  {
    subset_df <- arp_df[arp_df$title == pop, -c(1:2)]
    dna_object <- as.DNAbin(t(subset_df))
    return(dna_object)
  }),unique(arp_df$title))
  nuc_div_results <- lapply(dnabin_list, nuc.div)
  return(nuc_div_results)
}

#' Compute haplotype diversity for each sample
#'
#' This function computes the haplotype diversity for each sample in a haploid dataset. It converts the input dataframe to `DNAbin` objects for each population, removing columns that contain only missing values. It then computes the haplotype diversity and returns the results.
#'
#' @param arp_df A dataframe where each row represents a sequence. The dataframe must contain at least the following columns:
#'   - `title` : sample names or identifiers.
#'   - columns two and following : genetic data in haploid format. Columns containing only missing values (NA) will be removed during processing.
#' @return A named list where each element contains the haplotype diversity for one sample. The names of the list elements correspond to the unique sample identifiers in the `title` column of the input dataframe.
#'
#' @details
#' The function performs the following steps:
#' 1. For each unique sample (identified by the `title` column), it subsets the dataframe to include only the rows corresponding to that sample.
#' 2. Removes columns that consist entirely of missing values (NA).
#' 3. Converts the subsetted dataframe to a `DNAbin` object using the `as.DNAbin` function.
#' 4. Calculates haplotype diversity for each `DNAbin` object using the `hap.div` function with the "Nei" method.
#'
#' @examples
#' hap_div_results <- haploid.hap.div(arp_df)
#' print(hap_div_results)
#'
#' @export
haploid.hap.div <- function(arp_df)
{
  # Create a list of DNAbin objects for each population, removing columns with only NAs
  dnabin_list <- setNames(lapply(unique(arp_df$title), function(pop)
  {
    subset_df <- arp_df[arp_df$title == pop, -1]
    na_columns <- sapply(subset_df, function(col) all(is.na(col)))
    subset_df <- subset_df[, !na_columns]
    dna_object <- as.DNAbin(t(subset_df))
    return(dna_object)
  }), unique(arp_df$title))
  hap_div_results <- lapply(dnabin_list, hap.div, method="Nei")
  return(hap_div_results)
}

#' Compute pairwise Fst for each sample
#'
#' This function computes the pairwise Fst values for each pair of samples in a haploid dataset. It uses the `pairwise.neifst` function to compute Fst values.
#'
#' @param arp_df A dataframe where each row represents a sequence. The data frame must contain at least the following columns:
#'   - `title` : sample names or identifiers.
#'   - columns two and following : genetic data in haploid format.
#' @return A matrix containing the pairwise Fst values for each pair of samples. The rows and columns of the matrix correspond to the sample names, and each cell represents the Fst value between the sample pair.
#'
#' @details
#' The function calls `pairwise.neifst` with the parameter `diploid=F` to compute the Fst values for haploid data.
#'
#' @examples
#' fst_results <- haploid.fst(arp_df)
#' print(fst_results)
#'
#' @export
haploid.fst <- function(arp_df)
{
  fst_results <- pairwise.neifst(arp_df,diploid=F)
  return(fst_results)
}
