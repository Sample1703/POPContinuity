#Function to extract the number of samples
extract.num.samples <- function(folder)
{
  sam_file <- list.files(folder,pattern="\\.sam") # search for .SAM files in the settings folder
  file_content <- readLines(file.path(folder,sam_file))[-1] # read the .SAM file and remove the first line (corresponding to the number of samples)
  unique_keywords <- unique(gsub("\\d+$","",unique(sub("\\t.*","",file_content)))) #remove the final digit and what comes after and select the unique values
  num_samples <- length(unique_keywords) # the number of unique keywords
  return(num_samples)
}

## Split the simulated FASTA files
split.fasta.file.indiv.sim <- function(fasta_file)
{
  #Read the lines from the input FASTA file
  fasta_lines <- readLines(fasta_file)
  #Extract unique individual ID from lines starting with '>'
  ids <- unique(sub("^>([^_]+_[^_]+)_.*", "\\1", fasta_lines[startsWith(fasta_lines, ">")]))
  populations <- setNames(vector("list",length(ids)),ids)
  # Function to process lines based on id
  process.lines <- function(id)
  {
    result <- c()
    i <- 1
    while (i <= length(fasta_lines))
    {
      if (startsWith(fasta_lines[i],paste0(">",id)))
      {
        result <- c(result,fasta_lines[i],fasta_lines[i+1])
        i <- i+2  # Move to the next pair
      }
      else
      {
        i <- i+1  # Move to the next line if no match
      }
    }
    return(result)
  }

  # Apply the  function to each id
  populations <- lapply(ids,process.lines)
  names(populations) <- ids

  # Write the lines into different files for each population
  sapply(ids, function(id)
  {
    file_name <- paste0(id,".fasta")
    writeLines(populations[[id]],file_name)
  })
}

#Regrouper les individus par population
split.fasta.file <- function(fasta_file)
{
  #Read the lines from the input FASTA file
  fasta_lines <- readLines(fasta_file)
  #Extract unique population ID from lines starting with '>'
  ids <- unique(sub("^>([^_]+)_.*", "\\1",fasta_lines[startsWith(fasta_lines, ">")]))
  # Initialize a list for each population
  populations <- setNames(vector("list",length(ids)),ids)
  # Function to process lines based on id
  process.lines <- function(id)
  {
    result <- c()
    i <- 1
    while (i <= length(fasta_lines))
    {
      if (startsWith(fasta_lines[i],paste0(">",id)))
      {
        result <- c(result,fasta_lines[i],fasta_lines[i+1])
        i <- i+2  # Move to the next pair
      }
      else
      {
        i <- i+1  # Move to the next line if no match
      }
    }
    return(result)
  }

  # Apply the  function to each id
  populations <- lapply(ids,process.lines)
  names(populations) <- ids

  # Write the lines into different files for each population
  sapply(ids, function(id)
  {
    file_name <- paste0("population_",id, ".fasta")
    writeLines(populations[[id]],file_name)
  })
}

#' Create a Haplotype DataFrame for Vanessa Atalanta Application
#'
#' This function generates a haplotype dataframe from a list of FASTA files specifically for the Vanessa Atalanta application. It reads genetic data from the provided FASTA files, extracts SNP information, and organizes it into a dataframe.
#'
#' @param fasta_files_list A character vector containing file paths to the FASTA files with genetic sequences. Each file represents genetic data for one individual.
#'
#' @return A dataframe containing haplotype information. The dataframe includes:
#' \itemize{
#'   \item \strong{Individual}: The ID of the individual, derived from the file name.
#'   \item \strong{Population}: The population group to which the individual belongs, labeled based on specific years ("1910", "1950", "1960", "1970", and "2022").
#'   \item \strong{1 to 301}: SNP positions for each individual, representing haplotype information. Each SNP is captured from two sets of positions (1-301 and 302-602).
#' }
#'
#' @details
#' The function processes each FASTA file to extract SNPs for each individual. It creates a dataframe with columns for individual names, population labels, and SNP haplotypes. The populations are assigned labels as follows:
#' \itemize{
#'   \item \strong{1910}: First 10 individuals
#'   \item \strong{1950}: Next 9 individuals
#'   \item \strong{1960}: Next 10 individuals
#'   \item \strong{1970}: Next 7 individuals
#'   \item \strong{2022}: Last 7 individuals
#' }
#' The SNPs are combined from two different sets of positions into a single haplotype.
#'
#' @examples
#' \dontrun{
#' # Assuming fasta_files_sim is a vector of file paths to the FASTA files
#' vatalanta_df <- build.df.vatalanta(fasta_files_sim)
#' }
#'
#' @export
build.df.vatalanta <- function(fasta_files_list)
{
  #Create an empty dataframe
  haplotype_df_sim <- data.frame(matrix(ncol=302,nrow=length(fasta_files_list))) #create an empty dataframe
  colnames(haplotype_df_sim) <- c("Individual",1:301) #update the column names
  haplotype_df_sim$Individual <- gsub(".fasta","",basename(fasta_files_list))
  haplotype_df_sim$Population <- c(1:43)
  haplotype_df_sim$Population[1:10] <- "1910"
  haplotype_df_sim$Population[11:19] <- "1950"
  haplotype_df_sim$Population[20:29] <- "1960"
  haplotype_df_sim$Population[30:36] <- "1970"
  haplotype_df_sim$Population[37:43] <- "2022"
  haplotype_df_sim <- haplotype_df_sim[,c(1,303,2:302)] #reorder the column
  haplotype_df_sim <- haplotype_df_sim[c(1,3:10,2,11:20,22:29,21,30:43),] #reorder the rows

  snp_list_sim <- list() #create an empty list to store SNPs from all files
  for (fasta_file in fasta_files_list)
  {
    sequences <- read.fasta(fasta_file) #read sequences from the FASTA file
    indiv_snp_list_sim <- list() #create an empty list to store SNPs for the current file
    for (i in 1:length(sequences))
    {
      snp <- sequences[[i]][1] #retrieve the 51th position
      indiv_snp_list_sim[[i]] <- snp
      names(indiv_snp_list_sim)[i] <- positions[i] #add the position as the element name in the loop
    }
    snp_list_sim[[fasta_file]] <- indiv_snp_list_sim
  }
  names(snp_list_sim) <- gsub(".fasta", "",basename(names(snp_list_sim))) #update the name of the entry in the list

  split_snp_list_sim <- list()
  for (name in names(snp_list_sim))
  {
    #Split each element of snp_list_sim in two
    snp1 <- snp_list_sim[[name]][1:301]
    snp2 <- snp_list_sim[[name]][302:602]

    # Create a data frame to store the information
    df <- data.frame(
      position = names(snp1),
      haplotype = paste0(unlist(snp1),unlist(snp2)),
      stringsAsFactors = F)

    # Add the data frame to the result list
    split_snp_list_sim[[name]] <- df
  }
  individual_names <- names(split_snp_list_sim)
  for (i in seq_along(split_snp_list_sim))
  {
    individual_name <- individual_names[i]
    df <- split_snp_list_sim[[individual_name]]
    haplotype_df_sim[haplotype_df_sim$Individual==individual_name,3:303] <- unlist(df$haplotype)
  }
  return(haplotype_df_sim)
}

