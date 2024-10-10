#' Extract the mutation rate from a genetic file
#'
#' This function extracts the mutation rate from a genetic file. It reads the file located in the specified folder and retrieves the mutation rate.
#'
#' @param folder A character string specifying folder where the genetic file is located.
#' @return A numeric value representing the mutation rate.
#'
#' @details
#' The function uses the `extract.gen.file` function to get the name of the file containing the mutation rate. It then reads the file and extracts the mutation rate from the fourth element of the fourth line, assuming the mutation rate is positioned in this location in the file.
#'
#' @examples
#' extract.muta("example_folder")
#'
#' @export
extract.muta <- function(folder)
{
  genfile <- extract.gen.file(folder)
  muta <- strsplit(readLines(file.path(folder,genfile))," ")[[4]][4]
  return(muta)
}

#' Extract genetic files from a folder
#'
#' This function retrieves the names of genetic files with the `.par` extension from a specified folder.
#'
#' @param folder A character string specifying the folder where the `.par` files are located.
#' @return A character vector containing the names of the `.par` files found in the specified folder.
#'
#' @examples
#' extract.gen.file("example_folder")
#'
#' @export
extract.gen.file <- function(folder)
{
  genfiles <- list.files(paste0(folder,"/"),pattern=".par") # extract the vegK files where the value of K must be replaced
  return(genfiles)
}

#' Create a prior distribution for mutation rate
#'
#' This function generates a sequence of values representing a prior distribution for the mutation rate. The sequence ranges from 0.000002 to 0.000003 with a specified increment.
#'
#' @return A numeric vector containing the sequence of mutation rate values.
#'
#' @details
#' The function creates a sequence of values starting from 0.000002 and ending at 0.000003, with an increment of 0.00000001. This range is used as a prior distribution for the mutation rate in simulations.
#'
#' @examples
#' prior <- create.prior.muta()
#'
#' @export
create.prior.muta <- function()
{
  rangemuta <- seq(from = 0.000002, to = 0.000003, by = 0.00000001)
  return(rangemuta)
}

#' Extract a random mutation rate from a prior distribution
#'
#' This function samples a single random value from a given range of possible mutation rates.
#'
#' @param rangemuta A numeric vector representing the prior distribution of mutation rates.
#' @return A numeric value representing the randomly sampled mutation rate.
#'
#' @details
#' The function uses the `sample` function to randomly select one value from the provided numeric vector. This is used for simulations that require a random mutation rate drawn from a predefined prior distribution.
#'
#' @examples
#' prior_muta <- create.prior.muta()
#' random_muta <- extract.random.muta(prior_muta)
#'
#' @export
extract.random.muta <- function(rangemuta)
{
  randommuta <- sample(rangemuta, 1) # sample a single random value from the entire sequence
  return(randommuta)
}

#' Generate a list of random mutation rates
#'
#' This function creates a list of random mutation rate values by sampling from a given prior distribution. The length of the list is specified by the user.
#'
#' @param rangemuta A numeric vector representing the prior distribution of mutation rates.
#' @param nmb_simul An integer specifying the number of random mutation rates to generate, equal to the number of simulations to perform.
#' @return A numeric vector containing randomly sampled mutation rates from the provided sequence. The length of the vector equals the number of simulations specified.
#'
#' @details
#' The function uses the `extract.random.muta` function to sample mutation rates from the `rangemuta` vector. It repeats this sampling process `nmb_simul` times to create a list of random mutation rates.
#'
#' @examples
#' prior_muta <- create.prior.muta()
#' mutation_rates <- generate.muta.list(prior_muta,10)
#'
#' @export
generate.muta.list <- function(rangemuta,nmb_simul)
{
  muta_values <- replicate(nmb_simul,extract.random.muta(rangemuta))
  return(muta_values)
}

#' Replace mutation rate in genetic files
#'
#' This function updates the mutation rate in genetic file by overwriting the existing mutation rate with new value.
#'
#' @param folder A character string specifying the folder containing the genetic file.
#' @param muta_values A numeric vector containing the new mutation rate to be used in the genetic file.
#'
#' @return None. This function modifies the content of the genetic parameter files directly.
#'
#' @details
#' The function first extracts the genetic file name from the specified folder using the `extract.gen.file` function. It then reads the content of the file, replaces the existing mutation rate with as value from the `muta_values` vector, and writes the updated content back to the file.
#'
#' @examples
#' muta_values <- generate.muta.list(rangemuta,5)
#' replace.gen.file("example_folder",muta_values)
#'
#' @export
replace.gen.file <- function(folder,muta_values)
{
  genfile <- extract.gen.file(folder)
  content <- readLines(file.path(folder,genfile))
  for (i in seq_along(muta_values[i]))
  {
    newmuta <- muta_values[i]
    oldmuta <- strsplit(readLines(file.path(folder,genfile))," ")[[4]][4]
    content <- gsub(oldmuta,newmuta,content)
  }
  writeLines(content,file.path(folder,genfile))
}
