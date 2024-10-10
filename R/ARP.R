#' Extract content from ARP files
#'
#' This function searches for `.arp` files in the `GeneticsOutput` subdirectory of the specified folder and reads the content of the `.arp` file found.
#'
#' @param folder A character string specifying the name of the working directory folder containing the `GeneticsOutput` subdirectory.
#' @return A character vector containing the lines of the `.arp` file.
#' @examples
#' \dontrun{
#' content <- extract.ARP.files("example_folder")
#' }
#' @export
extract.ARP.files <- function(folder)
{
  arp_file <- list.files(paste0(folder,"/GeneticsOutput/"),pattern=".arp")
  content_arp <- readLines(paste0(folder,"/GeneticsOutput/",arp_file))
  return(content_arp)
}

#' Split a FASTA file by population IDs
#'
#' This function reads a FASTA file and splits it into separate FASTA files for each unique population based on identifiers in the sequence headers.
#'
#' @param fasta_file The path to the input FASTA file. The FASTA file should have sequence headers in the format ">ID_x", where "ID" is the unique identifier for the population.
#'
#' @return None. The function outputs multiple FASTA files, each corresponding to a unique population. These files are named "population_ID.fasta", where "ID" is the population identifier.
#'
#' @details
#' The function processes each line of the input FASTA file, extracts the population identifier from the sequence headers, and groups sequences by these identifiers. Each individual's sequences are written to a new FASTA file.
#'
#' @examples
#' split.fasta.file("example.fasta")
#'
#' @export
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
      if (startsWith(fasta_lines[i],paste0(">",id))) #if the line starts with ">"
      {
        result <- c(result,fasta_lines[i],fasta_lines[i+1]) #save the line and the next
        i <- i+2 #move to the next pair of lines
      }
      else
      {
        i <- i+1 #else, move to the next line
      }
    }
    return(result)
  }

  #Apply the  function to each id
  populations <- lapply(ids,process.lines)
  names(populations) <- ids

  #Write the lines into different files for each population
  sapply(ids, function(id)
  {
    file_name <- paste0("population_",id, ".fasta")
    writeLines(categories[[id]],file_name)
  })
}
