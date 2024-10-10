#' Build an empty dataframe for statistical analysis
#'
#' This function creates an empty data frame with the appropriate number of columns and rows based on the number of simulations and samples.
#'
#' @param folder A character string specifying the working directory containing the data.
#' @param nmb_simul An integer indicating the number of simulations.
#' @return A data frame with the specified dimensions, ready to contain the output of the statistical analysis.
#' @examples
#' \dontrun{
#' folder <- "example_folder"
#' nmb_simul <- 100
#' stat_df <- build.stat.df(folder,nmb_simul)
#' }
#' @export
build.stat.df <- function(folder,nmb_simul)
{
  num_samples <- extract.num.samples(folder) #retrieve the number of samples
  stat_df <- data.frame(matrix(ncol=(4+(num_samples*2)+num_samples*(num_samples-1)/2),nrow=nmb_simul+1)) #create a dataframe with the desired dimensions
  return(stat_df)
}

#' Save a dataframe to a .TXT file
#'
#' This function saves a data frame into a .TXT file named Stats_[folder].txt in the specified folder.
#'
#' @param folder A character string specifying the name of the working directory to be included in the output file name.
#' @param df A data frame to be saved to the .TXT file.
#' @return None. The function saves the dataframe to a .TXT file.
#' @examples
#' \dontrun{
#' folder <- "example_folder"
#' df <- data.frame(a = 1:5, b = letters[1:5])
#' df.to.file(folder,df)
#' }
#' @export
df.to.file <- function(folder,df)
{
  filename <- paste0("Stats_",folder,".txt") #save the name of the .txt file
  write.table(df,file=filename,sep="\t",row.names=F,col.names=T,quote=F,append=T)
}

#' Fill the content of the statistics dataframe
#'
#' This function computes nucleotide diversity and Fst statistics from the given ARP data frame and fills them into the provided statistics data frame.
#'
#' @param stat_df An empty dataframe to be filled with statistics results.
#' @param arp_df A dataframe containing genetic data from which statistics are computed.
#' @return The dataframe filled with computed statistics.
#' @examples
#' \dontrun{
#' stat_df <- data.frame(matrix(ncol = 3, nrow = 2)) # Example empty data frame
#' arp_df <- data.frame(title = c("sample1","sample2"), locus1 = c("A","T"), locus2 = c("C","G"), locus3 = c("A","C")) # Example ARP data frame
#' filled_df <- fill.stat.df(stat_df,arp_df)
#' }
#' @export
fill.stat.df <- function(stat_df,arp_df)
{
  #Compute the statistics
  nuc_div_results <- haploid.nuc.div(arp_df)
  #print(nuc_div_results)
  fst_results <- haploid.fst(arp_df)
  #print(fst_results)

  #Retrieve information on the samples
  num_samples <- length(unique(arp_df$title))
  sample_names <- unique(arp_df$title)

  #Add colnames to stat_df
  sample_combinations <- combn(sample_names,2,function(x) paste0("Fst",paste(x,collapse="/")),simplify=T) #create a combination of : "Fst",sample name,"/",sample name
  colnames(stat_df) <- c("SimID","K1","K2","MutationRate",paste0("nucdiv",sample_names),sample_combinations) #rename the columns

  #Fill the table with the nuc.div information
  stat_df[i,5:(as.numeric(num_samples)+3)] <- unlist(round(as.numeric(nuc_div_results),7)) #fill the columns 3 > num_samples+1 with the nucdiv variance values

  #Fill the table with the Fst information
  for (fst_column_name in colnames(stat_df)[grepl("^Fst",colnames(stat_df))]) #select the column names containing "Fst"
  {
    pair <- unlist(strsplit(sub("^Fst\\s*","",fst_column_name),"/")) #separate the pair at the "/" level (i.e. modern/neol1 -> modern + neol1)
    sample1 <- pair[1] #sample1 is the one before the "/"
    sample2 <- pair[2] #sample2 is the one after the "/"
    fst_value <- fst_results[sample1,sample2] #take the Fst value for the pair of interest
    stat_df[i,fst_column_name] <- round(as.numeric(fst_value),7) #put the value in the right line and right column
  }
  return(stat_df)
}
