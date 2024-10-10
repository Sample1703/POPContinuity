#' Create subsets of observed data for the ABC Analysis
#'
#' This function generates subsets of observed data from a dataframe to use in Approximate Bayesian Computation (ABC) analysis. The subsets are created based on pairwise combinations of haploid diversity columns and the corresponding Fst values.
#'
#' @param df A dataframe containing observed data with columns for haploid diversity and Fst values. The columns related to haploid diversity should be prefixed with "hapdiv" and the columns related to Fst values should be prefixed with "Fst".
#' @return A list of dataframes, each representing a subset of the observed data. Each subset dataframe contains columns for K1, the haploid diversity value for each population in the pair, and the Fst value between the population pair.
#'
#' @details
#' The function extracts column names related to haploid diversity from the dataframe and creates pairwise combinations of these columns. For each pair, it constructs a subset dataframe that includes the columns "K1", the two haploid diversity columns, and the corresponding Fst value column. Each subset is added to a list, which is then returned.
#'
#' @examples
#' subsets <- create.subsets.obs(df)
#'
#' @export
create.subsets.obs <- function(df)
{
  column_names <- colnames(df)
  intradiv_columns <- column_names[grep("^intra",column_names)]
  subsets <- list()
  for (i in 1:(length(intradiv_columns)-1))
  {
    for (j in (i+1):length(intradiv_columns))
    {
      fst_column <- paste0("Fst", sub("^intra","",intradiv_columns[i]),".",sub("^intra","",intradiv_columns[j]),sep="")
      subset_columns <- c("K1",intradiv_columns[i],intradiv_columns[j],fst_column)
      #print(subset_columns)
      subset_data <- df[,subset_columns,drop=F]
      subsets[[length(subsets) + 1]] <- subset_data
    }
  }
  return(subsets)
}

#' Create subsets of the main dataframe
#'
#' This function generates subsets from the main dataframe for each population pair. Each subset includes columns for simulation ID, carrying capacity value, nucleotide diversity, haploid diversity, and Fst values.
#'
#' @param df A dataframe containing columns related to nucleotide diversity, haploid diversity, and Fst values. The columns related to nucleotide diversity should be prefixed with "nucdiv", and those related to haploid diversity with "hapdiv". The dataframe should also include columns for "SimID" and "K1".
#' @return A list of dataframes, each representing a subset of the main dataframe. Each subset dataframe contains columns for simulation ID ("SimID"), carrying capacity value ("K1"), the haploid diversity value for each population in the pair, the nucleotide diversity value for each population in the pair, and the Fst value between the population pair.
#'
#' @details
#' The function extracts column names related to nucleotide diversity and haploid diversity from the dataframe. It then creates pairwise combinations of the nucleotide diversity columns. For each pair, it constructs a subset dataframe that includes the columns "SimID", "K1", the two nucleotide diversity columns, the two haploid diversity columns, and the Fst value column. Each subset is added to a list, which is then returned.
#'
#' The naming convention for the Fst columns is derived from the nucleotide diversity columns. For example, if the nucleotide diversity columns are "nucdiv1" and "nucdiv2", the corresponding Fst column will be named "Fst1/2".
#'
#' @examples
#' subsets <- create.subsets(df)
#'
#' @export
create.subsets <- function(df)
{
  column_names <- colnames(df)
  nucdiv_columns <- column_names[grep("^nucdiv",column_names)]
  hapdiv_columns <- column_names[grep("^hapdiv",column_names)]
  subsets <- list()
# Loop through intra columns
for (i in 1:(length(intra_columns)-1))
{
  for (j in (i+1):length(intra_columns))
  {
    fst_column <- paste0("Fst", sub("^intra","",intra_columns[i]),"/",sub("^intra","",intra_columns[j]),sep="")
    subset_columns <- c("SimID","K1",intra_columns[i],intra_columns[j],fst_column)
    subset_data <- stat_df_v2[,subset_columns,drop=F]
    subsets[[length(subsets) + 1]] <- subset_data
  }
}
  return(subsets)
}

#' Create subsets of the main dataframe (version V. Atalanta)
#'
#' This function generates subsets from the main dataframe for each population pair. Each subset includes columns for simulation ID, carrying capacity value, nucleotide diversity, haploid diversity, and Fst values.
#'
#' @param df A dataframe containing columns related to nucleotide diversity, haploid diversity, and Fst values. The columns related to nucleotide diversity should be prefixed with "nucdiv", and those related to haploid diversity with "hapdiv". The dataframe should also include columns for "SimID" and "K1".
#' @return A list of dataframes, each representing a subset of the main dataframe. Each subset dataframe contains columns for simulation ID ("SimID"), carrying capacity value ("K1"), the haploid diversity value for each population in the pair, the nucleotide diversity value for each population in the pair, and the Fst value between the population pair.
#'
#' @details
#' The function extracts column names related to nucleotide diversity and haploid diversity from the dataframe. It then creates pairwise combinations of the nucleotide diversity columns. For each pair, it constructs a subset dataframe that includes the columns "SimID", "K1", the two nucleotide diversity columns, the two haploid diversity columns, and the Fst value column. Each subset is added to a list, which is then returned.
#'
#' The naming convention for the Fst columns is derived from the nucleotide diversity columns. For example, if the nucleotide diversity columns are "nucdiv1" and "nucdiv2", the corresponding Fst column will be named "Fst1/2".
#'
#' @examples
#' subsets <- create.subsets(df)
#'
#' @export
create.subsets.alt <- function(df) {
  column_names <- colnames(df)

  # Identify intra and Fst columns
  intra_columns <- column_names[grep("^intra", column_names)]

  subsets <- list()

  # Loop through intra columns
  for (i in 1:(length(intra_columns)-1))
  {
    for (j in (i+1):length(intra_columns))
    {
      fst_column <- paste0("Fst", sub("^intra","",intra_columns[i]),"/",sub("^intra","",intra_columns[j]),sep="")
      subset_columns <- c("SimID","K1",intra_columns[i],intra_columns[j],fst_column)
      subset_data <- df[,subset_columns,drop=F]
      subsets[[length(subsets) + 1]] <- subset_data
    }
  }
  return(subsets)
}

#' Save a data subset to a text file
#'
#' This function saves a subset of data to a tab-separated text file. The file is named based on the provided index to distinguishing between the different subsets.
#'
#' @param subset A dataframe representing a subset of the main dataframe. This dataframe will be written to a text file.
#' @param index A list used to name the output file. The file will be named `subset_<index>.txt`, where `<index>` is replaced by the provided index value from the list.
#' @return None. The function writes the subset dataframe to a text file.
#'
#' @examples
#' index <- as.list(1:10)
#' subset.to.file(subset,index)
#'
#' @export
subset.to.file <- function(subset,index)
{
  filename <- paste0("subset_",index,".txt")
  write.table(subset,file=filename,sep="\t",quote=F,row.names=F)
}
