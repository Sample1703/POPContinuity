#' Calculate expected heterozygosity (Hs) for a genind Object
#'
#' This function calculates the expected heterozygosity (Hs) for a given genind object.
#'
#' @param genind A `genind` object.
#'
#' @return A numeric vector. The function returns the expected heterozygosity (Hs) values for the populations in the `genind` object.
#'
#' @details
#' Expected heterozygosity (Hs) is a measure of the intra-population diversity. The function uses the `Hs` function from the `adegenet` package to compute these values.
#'
#' @examples
#' \dontrun{
#' Hs_results <- Hs.div(genind)
#' }
#'
#' @importFrom adegenet Hs
#' @export
Hs.div <- function(genind)
{
  Hs_results <- Hs(genind)
  return(Hs_results)
}

#' Calculate the Nei Fst for a genind Object
#'
#' This function calculates the Fst values for a given genind object using Nei's method.
#'
#' @param genind A `genind` object.
#'
#' @return A matrix. The function returns a matrix of pairwise Fst values between populations.
#'
#' @examples
#' \dontrun{
#' fst_results <- diploid.fst(genind)
#' }
#'
#' @importFrom hierfstat genet.dist
#' @export
diploid.fst <- function(genind)
{
  fst_results <- genet.dist(genind,method="Nei87")
  return(fst_results)
}
