#' Launch a SPLATCHE3 simulation
#'
#' This function launches a SPLATCHE3 simulation using a job script. It ensures that the SPLATCHE3 executable has the correct permissions and then submits the job using the `sbatch` command.
#'
#' @param scriptname A character string specifying the name of the script to be executed to launch a simulation.
#' @param settingsfile A character string specifying the path to the settings file used by SPLATCHE3. This file should contain the parameters and settings required for the simulation.
#' @return A SLURM file and a log file. The SLURM file include information about the job status and the log file contains details on whether the simulation was successful.
#'
#' @details
#' The function first makes the SPLATCHE3 executable (`SPLATCHE3-Linux-64b`) runnable by changing its permissions. It then submits a job to the scheduler using `sbatch`, specifying the `fn-mat` partition.
#'
#' @examples
#' launch.simulation("script.sh","settings.txt")
#'
#' @export
launch.simulation <- function(scriptname,settingsfile) # function to launch 1 simulation
{
  system("chmod +x SPLATCHE3-Linux-64b") # make the file executable
  system(paste("sbatch -p fn-mat",scriptname,settingsfile)) # run the script
}

#' Check for the existence of the `GeneticsOutput` folder
#'
#' This function checks whether the `GeneticsOutput` folder exists within the working directory. If the folder is not found, the function waits for 1 second and then recursively checks again until the folder is found.
#'
#' @param folder A character string specifying the working directory inside which the `GeneticsOutput` folder should be located.
#' @return A logical value. Returns `TRUE` if the `GeneticsOutput` folder exists; otherwise, it continues checking until the folder is found.
#'
#' @examples
#' check.for.file("example_folder")
#'
#' @export
check.for.file <- function(folder)
{
  if (dir.exists(paste0(folder,"/GeneticsOutput")))
  {
    return(TRUE)
  }
  else
  {
    Sys.sleep(1)
    return(check.for.file(folder))
  }
}
