#' Extract K and m values from settings file
#'
#' This function extracts the values of carrying capacity (K) from vegetation files and migration rate (m) from a settings text file within the specified working directory
#'
#' @param folder A character string specifying the folder containing the settings and vegetation files.
#' @return A list containing the migration rate (m) as a numeric value and the carrying capacity values (K) as a list of integers.
#' @examples
#' \dontrun{
#' folder <- "example_folder"
#' Km_values <- extract.Km(folder)
#' print(Km_values)
#' }
#' @export
extract.Km <- function(folder)
{
  vegKfiles <- extract.vegK.files(folder) # retrieve the names of the vegK files
  valuesK <- list()
  for (file in vegKfiles)
  {
    valueK <- strsplit(readLines(file.path(folder, file))[2],"\t")[[1]][2] # extract the value of K from the file
    valueK <- as.integer(valueK) #convert it as an integer
    valuesK <- c(valuesK, valueK) # append the value of K to the list
  }
  valuem <- na.omit(str_match(readLines(paste0(folder,".txt")),"MigrationRate=([0-9.]+)"))[2] # extract the value of m
  valuesKm <- list(as.numeric(valuem), valuesK) # save the values of m and K in a list
  return(valuesKm)
}

#' Create prior distributions for carrying capacity (K) values
#'
#' This function generates prior distributions for carrying capacity values (K) based on the values provided. Each prior is a sequence of values ranging from K/10 to K*10.
#'
#' @param valueKm A list containing migration rate (m) and a list of carrying capacities (K) as returned by `extract.Km()`.
#' @return A list of vectors, each representing a range of prior values for the corresponding K.
#' @examples
#' \dontrun{
#' valueKm <- list(0.01,list(50,100))
#' priors <- create.prior(valueKm)
#' print(priors)
#' }
#' @export
create.prior <- function(valueKm)
{
  rangeK <- list() # create a list that will contain the priors
  for (valueK in valueKm[[2]]) # loop for each value of K retrieved using extract.Km()
  {
    range <- seq(round(valueK / 10,0),round(valueK * 10,0)) # create a vector of values ranging from K/10 to K*10
    rangeK[[length(rangeK)+1]] <- range # append the new prior at the end of the list of priors
  }
  return(rangeK)
}

#' Extract a random carrying capacity value (K) from a prior
#'
#' This function extracts a random value of carrying capacity (K) from each range provided.
#'
#' @param rangeK A list of vectors, each representing a range of prior values for K.
#' @return A list of single random values of K, one from each prior range.
#' @examples
#' \dontrun{
#' rangeK <- list(seq(5,500,1), seq(10,1000,1))
#' randomK <- extract.random.K(rangeK)
#' print(randomK)
#' }
#' @export
extract.random.K <- function(rangeK)
{
  randomK <- lapply(rangeK, function(x) sample(x,1)) # sample a single random value of K per prior
  return(randomK)
}

#' Generate a list of random carrying capacity (K) values
#'
#' This function generates a list of random carrying capacity (K) values from provided ranges. It replicates the process of selecting a random K value from each range a specified number of times.
#'
#' @param rangeK A list where each element is a vector of possible values for K. Each vector represents a range from which random K values will be drawn.
#' @param nmb_simul An integer specifying the number of random K values to generate for each range (one per simulation).
#' @return A list of length `nmb_simul`, where each element is a vector of randomly sampled K values within the ranges specified in `rangeK`.
#'
#' @examples
#' ranges <- list(seq(5,500,1), seq(10,1000,1))
#' K_list <- generate.K.list(ranges,10)
#' print(K_list)
#'
#' @export
generate.K.list <- function(rangeK,nmb_simul)
{
  K_values <- replicate(nmb_simul,extract.random.K(rangeK))
  return(K_values)
}

#' Extract names of vegetation (VegK) files from a folder
#'
#' This function reads the `dynamic_K.txt` file from a specified folder and extracts the names of the vegetation (VegK) files listed within. It assumes that `dynamic_K.txt` contains a list of file paths where the first line indicates the number of files, and the following lines list their paths.
#'
#' @param folder A character string specifying the path to the folder containing the `dynamic_K.txt` file.
#' @return A character vector of filenames for VegK files extracted from the `dynamic_K.txt` file. The path is removed from the filenames, keeping only the filenames with their `.txt` extension.
#'
#' @examples
#' vegK_files <- extract.vegK.files("example_folder")
#' print(vegK_files)
#'
#' @export
extract.vegK.files <- function(folder)
{
  dynamic_K_content <- readLines(file.path(folder, "dynamic_K.txt")) # read the content of the dynamic_K file
  vegKfiles <- dynamic_K_content[-1] # remove the first line containing the number of files
  vegKfiles <- gsub(".*/(\\S+\\.txt)\\s.*", "\\1", vegKfiles) # retrieve the names of the vegK files
  return(vegKfiles)
}

#' Replace the value of carrying capacity (K) in specified files
#'
#' This function replaces the value of carrying capacity (K) in the files listed in `dynamic_K.txt` within a specified folder. The function reads the content of each file, replaces occurrences of K with a new value, and then overwrites the updated content in the original files.
#'
#' @param folder A character string specifying the folder containing the `dynamic_K.txt` file and the vegetation (VegK) files where K should be replaced.
#' @param Kall A list of new K values to be used for replacement.
#' @details
#' The function reads the file paths from `dynamic_K.txt` in the given folder. It then updates the value of K in each file, where K is replaced in fields surrounded by tab characters.
#' @return None. This function modifies the contents of the VegK files.
#'
#' @examples
#' folder <- "example_folder"
#' new_K_values <- list(10,20,30)
#' replace.K.file(folder,new_K_values)
#'
#' @export
replace.K.file <- function(folder,Kall)
{
  vegKfiles <- extract.vegK.files(folder) # extract the files where the value of K must be replaced
  content <- readLines(file.path(folder,vegKfiles)) # read the content of the file
  Knew <- Kall[[i]]
  content <- gsub("(?<=\\t)\\d+(?=\\t)",Knew,content,perl=T) # look for one or more digits (in between 2 tabulations) and replace them with the new value of K
  writeLines(content, file.path(folder,vegKfiles)) # write the updated file (original file is replaced)
}


replace.K.file.alt <- function(folder,Kall)
{
  vegKfiles <- extract.vegK.files(folder) # extract the files where the value of K must be replaced
  content <- readLines(file.path(folder,vegKfiles)) # read the content of the file
  Knew <- Kall
  content <- gsub("(?<=\\t)\\d+(?=\\t)",Knew,content,perl=T) # look for one or more digits (in between 2 tabulations) and replace them with the new value of K
  writeLines(content, file.path(folder,vegKfiles)) # write the updated file (original file is replaced)
}
