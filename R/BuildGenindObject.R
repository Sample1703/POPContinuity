#' Build haploid genind object from ARP dataframe
#'
#' This function converts an ARP data frame into a haploid genind object using the `df2genind` function from the `adegenet` package.
#'
#' @param arp_df A data frame containing ARP data. The first column is assumed to be identifiers, and the following columns are haploid genetic data.
#' @return A genind object with haploid data.
#' @examples
#' \dontrun{
#' arp_df <- data.frame(
#'   ID = c("id1","id2"),
#'   locus1 = c("A","T"),
#'   locus2 = c("C","G")
#' )
#' genind_obj <- build.genind.haploid(arp_df)
#' }
#' @importFrom adegenet df2genind
#' @export
build.genind.haploid <- function(arp_df)
{
  genind_obj <- df2genind(arp_df[2:length(arp_df)],ploidy=1,ncode=1)
  return(haploid_genind_obj)
}

#' Build diploid genind object from ARP dataframe
#'
#' This function converts an ARP data frame into a diploid genind object using the `df2genind` function from the `adegenet` package.
#'
#' @param arp_df A data frame containing ARP data. The first column is assumed to be identifiers, and the following columns are diploid genetic data.
#' @return A genind object with diploid data.
#' @examples
#' \dontrun{
#' arp_df <- data.frame(
#'   ID = c("id1","id2"),
#'   locus1 = c("A/A","T/T"),
#'   locus2 = c("C/C","G/G")
#' )
#' genind_obj <- build.genind.diploid(arp_df)
#' }
#' @importFrom adegenet df2genind
#' @export
build.genind.diploid <- function(arp_file)
{
  genind_obj <- df2genind(arp_df[2:length(arp_df)],ploidy=2,ncode=2)
  return(diploid_genind_obj)
}
