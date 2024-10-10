build.abcdf <- function(stat_df)
{
  num_pairs <- length(grep("^Fst",names(stat_df),value=TRUE)) #retrieve the number of pairs
  pair_names <- unique(sub("^Fst","",names(stat_df[grep("^Fst",names(stat_df))]))) #retrieve unique pair names
  abc_df <- data.frame(matrix(ncol=6,nrow=num_pairs)) #create a dataframe with the desired dimensions
  colnames(abc_df) <- c("pair","intra1","intra2","Fst1/2","goodness of fit","p-value") #rename the columns
  abc_df[1] <- pair_names #fill the first column with the pair names
  return(abc_df)
}
