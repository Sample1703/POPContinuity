#' Build gene dataframe from ARP file
#'
#' This function processes an ARP file to extract sample names, sample sizes, and sequences, and then formats this information into a data frame suitable for analysis.
#'
#' @param arp_file A character vector containing the lines of an ARP file.
#' @return A data frame with the sample titles and their corresponding sequences split into individual characters.
#' @examples
#' \dontrun{
#' arp_file <- extract.ARP.files("example_folder")
#' gene_df <- build.genedf.kw1(arp_file)
#' }
#' @export
build.genedf.kw1 <- function(arp_file)
{
  #Retrieve the main information from the ARP file samples
  sample_names <- list((sub("^\\s*SampleName=\"(.*?)\"", "\\1",grep("^\\s*SampleName=\"(.+)\"",arp_file,value=TRUE)))) #sample names
  sample_sizes <- list(as.integer(gsub("^\\s*SampleSize=(\\d+).*", "\\1",grep("^\\s*SampleSize=\\d+",arp_file,value=TRUE)))) #sample sizes
  sample_info <- data.frame(sample_names,sample_sizes) #combine the information in a dataframe
  colnames(sample_info) <- c("Names", "Sizes") #rename the columns of the dataframe

  #Retrieve the sequences for each sample
  merged_list <- list() # create an empty list to store merged populations
  added_samples <- character(0) # create an empty set to store added sample names
  for (sample in unlist(sample_info$Names))
  {
    if (!(sample %in% added_samples)) # if the sample has not been added yet
    {
      merged_list[[sample]] <- select.lines(sample,arp_file) # select the line
      added_samples <- union(added_samples,sample) # add the sample to the set of added samples
    }
    else # if the sample has already been added
    {
      next # don't add it again to the set
    }
  }
  #Reformat the outputted list
  formatted_list <- format.list(merged_list)
  #Turn the list into a dataframe
  arp_df <- do.call(rbind, lapply(names(formatted_list), function(title)
  {
    keyword <- gsub("^(.*)_.*", "\\1", title)  # extract the keyword
    data.frame(title = keyword, sequence = unlist(formatted_list[[title]]),stringsAsFactors=F)
  }))
  #Split the sequences between each letter
  sequence_df <- as.data.frame(do.call(rbind,strsplit(arp_df$sequence, "")),stringsAsFactors=F)
  arp_df <- cbind(arp_df["title"],sequence_df)
  return(arp_df)
}

#' Build gene dataframe with modified titles from ARP file
#'
#' This function processes an ARP file to extract sample names, sample sizes, and sequences, and then formats this information into a data frame. Titles are modified by removing final digits.
#'
#' @param arp_file A character vector containing the lines of an ARP file.
#' @return A data frame with the sample titles and their corresponding sequences split into individual characters. Titles are modified by removing any final digits.
#' @examples
#' \dontrun{
#' arp_file <- extract.ARP.files("example_folder")
#' gene_df <- build.genedf.kw2(arp_file)
#' }
#' @export
build.genedf.kw2 <- function(arp_file)
{
  #Retrieve the main information from the ARP file samples
  sample_names <- list((sub("^\\s*SampleName=\"(.*?)\"", "\\1",grep("^\\s*SampleName=\"(.+)\"",arp_file,value=TRUE)))) #sample names
  sample_sizes <- list(as.integer(gsub("^\\s*SampleSize=(\\d+).*", "\\1",grep("^\\s*SampleSize=\\d+",arp_file,value=TRUE)))) #sample sizes
  sample_info <- data.frame(sample_names,sample_sizes) #combine the information in a dataframe
  colnames(sample_info) <- c("Names", "Sizes") #rename the columns of the dataframe

  #Retrieve the sequences for each sample
  merged_list <- list() # create an empty list to store merged populations
  added_samples <- character(0) # create an empty set to store added sample names
  for (sample in unlist(sample_info$Names))
  {
    if (!(sample %in% added_samples)) # if the sample has not been added yet
    {
      merged_list[[sample]] <- select.lines(sample,arp_file) # select the line
      added_samples <- union(added_samples,sample) # add the sample to the set of added samples
    }
    else # if the sample has already been added
    {
      next # don't add it again to the set
    }
  }
  #Reformat the outputted list
  formatted_list <- format.list(merged_list)
  #Turn the list into a dataframe
  arp_df <- do.call(rbind, lapply(names(formatted_list), function(title) {
    keyword <- gsub("\\d*$","", title)  # remove digits at the end of the title
    data.frame(title = keyword, sequence = unlist(formatted_list[[title]]), stringsAsFactors = FALSE)
  }))
  #Split the sequences between each letter
  sequence_df <- as.data.frame(do.call(rbind,strsplit(arp_df$sequence, "")),stringsAsFactors=F)
  arp_df <- cbind(arp_df["title"],sequence_df)
  #print(table(arp_df$title))
  return(arp_df)
}

#' Extract information from ARP file
#'
#' This function extracts and returns key information from an ARP file, including the title of the simulation, the number of samples, and details about each sample (names and sizes).
#'
#' @param arp_file A character vector containing the lines of an ARP file.
#' @return A named list containing the following elements:
#' \describe{
#'   \item{\code{Title}}{The title of the simulation extracted from the ARP file.}
#'   \item{\code{NumSamples}}{The number of samples specified in the ARP file.}
#'   \item{\code{Names}}{A list of sample names.}
#'   \item{\code{Size}}{A list of sample sizes.}
#' }
#' @examples
#' \dontrun{
#' arp_file <- extract.ARP.files("example_folder")
#' info <- extract.info(arp_file)
#' }
#' @export
extract.info <- function(arp_file)
{
  arp_title <- (sub("^\\s*Title=\"(.*?)\"", "\\1",grep("^\\s*Title=\"(.+)\"",arp_file,value=TRUE))) #title of the simulation
  num_samples <- as.list(as.integer(gsub("^\\s*NbSamples=(\\d+).*","\\1",grep("^\\s*NbSamples=\\d+",arp_file,value=TRUE)))) #number of samples
  sample_names <- as.list((sub("^\\s*SampleName=\"(.*?)\"", "\\1",grep("^\\s*SampleName=\"(.+)\"",arp_file,value=TRUE)))) #sample names
  sample_sizes <- as.list(as.integer(gsub("^\\s*SampleSize=(\\d+).*", "\\1",grep("^\\s*SampleSize=\\d+",arp_file,value=TRUE)))) #sample sizes
  info_list <- list("Title"=arp_title,"NumSamples"=num_samples,"Names"=sample_names,"Size"=sample_sizes)
  return(info_list)
}

#' Select lines for a specific sample from ARP file
#'
#' This function extracts lines corresponding to a specified sample from an ARP file. It identifies lines based on the sample name and captures those that start with a digit.
#'
#' @param sample_name A character string specifying the name of the sample to be selected.
#' @param arp_file A character vector containing the lines of an ARP file.
#' @return A character vector containing the lines from the ARP file that correspond to the specified sample and match the digit pattern.
#' @examples
#' \dontrun{
#' arp_file <- extract.ARP.files("example_folder")
#' sample_lines <- select.lines("modern",arp_file)
#' }
#' @export
select.lines <- function(sample_name,arp_file)
{
  current_sample <- NULL #pour conserver le nom du sample actuellement dans la boucle (i.e. "modern","neol",etc)
  capture_lines <- FALSE #pour déterminer si une ligne doit être capturée comme faisant partie du sample ou non (F=non,T=oui)
  selected_lines <- lapply(arp_file,function(line)
  {
    if (grepl("SampleName",line)) #si la ligne contient "SampleName"
    {
      current_sample <- gsub('"','',regmatches(line,regexpr('"([^"]*)"',line))[[1]]) #on extrait le nom du sample
      capture_lines <<- current_sample == sample_name #on compare le nom du sample actuel au nom du sample d'intérêt ; si identiques capture_lines est "TRUE"
    }
    if (capture_lines && grepl("^\\d+_\\d+",line)) #si capture_lines est "TRUE" et que la ligne commence par digit(s)_digit(s)
    {
      return(line)
    }
    else
    {
      return(NULL)
    }
  })
  selected_lines <- selected_lines[sapply(selected_lines,Negate(is.null))] #on sélectionne uniquement les lignes qui ne sont pas NULL
  return(selected_lines)
}

#' Format selected lines from ARP file
#'
#' This function formats a list of selected lines from an ARP file by splitting each line into parts based on tab delimiters and organizing them into a list. The names of the list elements correspond to the individual IDs extracted from the lines, and the content is the associated sequence.
#'
#' @param selected_lines_list A list of character vectors, where each character vector represents lines for a specific sample from an ARP file. Each line should be tab-separated.
#' @return A list of named lists, where each sublist contains the sequences for a particular sample. The names of the sublist elements are the individual IDs, and the values are the associated sequences.
#' @examples
#' \dontrun{
#' selected_lines <- list(
#'   modern = c("1_1\t1\tACGT", "1_2\t1\tTGCA"),
#'   ancient = c("2_1\t1\tCGTA", "2_2\t1\tGATC")
#' )
#' formatted_list <- format.list(selected_lines)
#' }
#' @export
format.list <- function(selected_lines_list)
{
  lapply(selected_lines_list, function(sublist) #on loop sur chaque sublist the selected_lines_list (i.e. "modern" et "lbk")
  {
    setNames(lapply(sublist, function(line) #on loop sur chaque ligne de la sublist (i.e. [[1]],[[2]],...)
    {
      parts <- strsplit(line,"\t")[[1]] #on split la ligne en séparant au niveau des tabulations
      element_name <- parts[1] #on récupère le 1er élément (ID de l'individu)
      element_content <- parts[3] #on récupère le 3ème élément (séquence associée)
      return(element_content)
    }),
    sapply(sublist, function(line) strsplit(line,"\t")[[1]][1]))
  })
}

#' Convert formatted list to dataframe
#'
#' This function converts a formatted list of sequences into a data frame. The list should be organized such that each population contains named lists of individual sequences. Each sequence is split into its constituent nucleotides, and a data frame is created where each row corresponds to a nucleotide for a specific individual and population.
#'
#' @param formatted_list A named list where each element represents a population and contains named lists of individual sequences. Each sequence is a character string of nucleotides.
#' @return A dataframe where each row corresponds to a nucleotide in the sequence for an individual within a population. The data frame includes columns for individual IDs, population names, and each nucleotide position.
#' @examples
#' \dontrun{
#' formatted_list <- list(
#'   modern = list(
#'     individual1 = "ACGT",
#'     individual2 = "TGCA"
#'   ),
#'   ancient = list(
#'     individual3 = "GATC",
#'     individual4 = "CTAG"
#'   )
#' )
#' df <- formatted.to.df(formatted_list)
#' }
#' @export
formatted.to.df <- function(formatted_list)
{
  do.call(rbind,lapply(names(formatted_list),function(population) #retrieve the population names in the list, apply the function to each population name and combine the results in a dataframe
  {
    do.call(rbind,lapply(names(formatted_list[[population]]),function(individual) #retrieve the names of the individuals and apply the function to each individual
    {data.frame(ID = individual, #create the column "ID" containing the name of the individuals
                Population = population, #create the column "Population" with the population associated with the individual
                t(matrix(strsplit(formatted_list[[population]][[individual]],"")[[1]]))) #split the sequence between each nucleotide
    }))
  }))
}

#' Save ARP dataframe and information to file
#'
#' This function saves an ARP data frame and an information list to a text file. The information list contains metadata about the ARP file, which is written as a title line, followed by the data frame written in tab-separated format.
#'
#' @param arp_df A dataframe containing ARP data.
#' @param info_list A list containing metadata about the ARP file. The first element of the list should be the title of the ARP file.
#' @return None.
#' @examples
#' \dontrun{
  #' arp_df <- data.frame(ID = c("id1","id2"), Population = c("pop1","pop2"), locus1 = c("A","T"), locus2 = c("C","G"))
#' info_list <- list("Sample Title",10,c("id1","id2"),c(2,2))
#' arpdf.to.file(arp_df,info_list)
#' }
#' @export
arpdf.to.file <- function(arp_df,info_list)
{
  write(paste0("Title line: \"",info_list[[1]],"\""),file="arp_data.txt") #write the name of the ARP file
  write.table(arp_df,file="ARP_summary",sep="\t",row.names=F,col.names=F,quote=F,append=T) #write the dataframe in the file
}
