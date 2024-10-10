#' Move files based on pattern
#'
#' This function moves files from their source folder to a destination folder if their filenames contain a specified pattern.
#' If the folder does not exist, it will be created.
#'
#' @param folder A character string specifying the path to the destination folder where the files should be moved.
#' @param pattern A character string representing the pattern to match in the filenames. Only files with names containing this pattern will be moved.
#' @return None. This function performs file operations, moving files that match the pattern to the specified folder.
#' @examples
#' \dontrun{
#' move.files("destination_folder","pattern")
#' }
#' @export
move.files <- function(folder,pattern)
{
  if (!file.exists(folder)) #if the folder does not exist
  {
    dir.create(folder) #create it
  }
  files <- list.files(pattern=pattern) #select the files whose names match the pattern
  for (file in files)
  {
    file.rename(file,file.path(folder,file)) #move the files in the folder
  }
}
