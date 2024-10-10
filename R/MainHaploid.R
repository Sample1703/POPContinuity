#' Main function to perform a population continuity test
#'
#' This function performs a series of computational steps, including setting up priors, running SPLATCHE3 simulations, computing genetic diversity statistics, performing ABC analysis, and conducting a continuity test.
#'
#' @param nmb_simul An integer corresponding to the number of simulations to perform.
#' @param folder A character string specifying the folder containing the necessary input files for the simulations and the statistical analysis.
#'
#' @return A dataframe with columns `Population`, `PValue`, and `Conclusion`. The `Population` column indicates the population pair, the `PValue` column provides the p-value from the comparison, and the `Conclusion` column states whether continuity is observed (`"Continuity"`) or not (`"No continuity"`).
#'
#' @details
#' The function follows these steps:
#' \enumerate{
#'   \item Retrieves the carrying capacity value provided by the user and creates a prior distribution.
#'   \item Creates a list of carrying capacity (K) values of length nmb_simul.
#'   \item Retrieves the mutation rate provided by the user and creates a prior distribution.
#'   \item Creates a list of mutation rate values of length nmb_simul.
#'   \item Initializes an empty data frame to store computed diversity statistics.
#'   \item For each K value, updates the values in the files, launches a SPLATCHE3 simulation, converts ARP files, and computes diversity statistics (nucleotide diversity, haplotype diversity, and Fst).
#'   \item Performs ABC analysis on each subset of simulated data and computes goodness-of-fit statistics with the associated p-values.
#'   \item Computes the mean of the best simulations for each statistic and pair.
#'   \item Conducts continuity test
#'   \item Cleans up log and subset files.
#'   \item Returns the results of the continuity test.
#' }
#'
#' @examples
#' folder <- "example_folder"
#' nmb_simul <- 1000
#' main.function(nmb_simul,folder)
#'
#' @export
main.haploid <- function(nmb_simul,folder)
{
  #[COMPUTATION AROUND K AND ITS PRIOR]

  #Retrieve the value of K provided by the user (Kuser)
  Km <- extract.Km(folder) #list of two elements containing the values of K and m
  Kuser <- Km[[2]] #value of Kuser
  #Create the prior and the list of nmb_simul K values
  Kprior <- create.Kprior(Km) #create a distribution of values ranging from Kuser/10 to Kuser*10
  Klist <- generate.K.list(Kprior,nmb_simul) #create a list of nmb_simul values of K taken randomly from the distribution
  Kall <- cbind(Kuser,Klist)  #combine the values of Kuser and randomly-drawn K in a list

  #[COMPUTATION AROUND MUTATION RATE AND ITS PRIOR]

  muta <- extract.muta(folder) #retrieve the value of m provided by the user
  mutaprior <- create.prior.muta() #create a mutation rate distribution
  mutalist <- generate.muta.list(mutaprior,nmb_simul) #create a list of nmb_simul values of m taken randomly from the distribution

  #[CREATION OF AN EMPTY DATAFRAME TO RETRIEVE THE COMPUTED STATISTICS]

  stat_df <- build.stat.df(folder,nmb_simul) #create the empty dataframe with the desired dimensions
  stat_df[1] <- c(1:(length(Kall)/2)) #fill the first column of the table with the simulation ID

  #[SIMULATION + REFORMAT ARP + STATISTICS COMPUTATION]

  #For each value of K is the list
  for (i in 1:(length(Kall)/2))
  {
    #UPDATE THE VALUES OF K AND MUTATION RATE IN TABLE AND FILES

    #Fill the table with the K and mutation rate values from the list previously created
    stat_df[i,2] <- Kall[1,i] #fill the 2nd column of the table with the K values
    stat_df[i,3] <- mutalist[i] #fill the 3rd column of the table with the mutation rate values

    #Identify the files where the value of K is present and update them with the one from the table
    replace.K.file(folder, Kall)
    #Identify the files where the value of mutation rate is present and update them with the one from the table
    replace.gen.file(folder,mutalist)

    #LAUNCH A SIMULATION

    launch.simulation("simulation.sh",paste0(folder,".txt")) #launch the simulation using the script simulation.sh and the setting file named as the settings folder

    #Wait until the folder "GeneticsOutput is created to continue the procedure
    if (check.for.file(folder))
    {
      arp_file <- extract.ARP.files(folder) #folder exists, extract the content of the ARP file it contains
    }
    else
    {
      cat("GeneticsOutput folder not created.")
    }

    #CONVERT ARP FILE

    arp_df <- build.genedf_newmerge(arp_file) #Convert the ARP file into a genetic dataframe

    #FILL THE STATISTICS TABLE

    #Compute the diversity statistics
    nuc_div_results <- haploid.nuc.div(arp_df) #nucleotide diversity
    hap_div_results <- haploid.hap.div(arp_df) #haplotype diversity
    fst_results <- haploid.fst(arp_df) #Fst

    #Retrieve information on the samples
    num_samples <- length(unique(arp_df$title)) #number of samples
    sample_names <- unique(arp_df$title) #names of the samples

    #Update the column names of the statistics table
    sample_combinations <- combn(sample_names,2,function(x) paste0("Fst",paste(x,collapse="/")),simplify=T) #create a combination of : "Fst",sample name,"/",sample name
    colnames(stat_df) <- c("SimID","K1","MutationRate",paste0("nucdiv",sample_names),paste0("hapdiv",sample_names),sample_combinations) #rename the columns

    #Fill the table with the nuc.div information
    stat_df[i,4:(as.numeric(num_samples)+3)] <- unlist(round(as.numeric(nuc_div_results),7)) #fill the columns 3 to num_samples+1 with the nucleotide diversity values

    #Fill the table with the hap.div information
    stat_df[i,(as.numeric(num_samples)+4):(2*(as.numeric(num_samples))+3)] <- unlist(round(as.numeric(hap_div_results),7)) #fill the columns num_samples+4 to 2*(num_samples+3) with the haplotype diversity values

    #Fill the table with the Fst information
    for (fst_column_name in colnames(stat_df)[grepl("^Fst",colnames(stat_df))]) #select the column names containing "Fst" in their name
    {
      pair <- unlist(strsplit(sub("^Fst\\s*","",fst_column_name),"/")) #separate the pair at the "/" level (i.e. modern/neol -> modern + neol)
      sample1 <- pair[1] #sample1 is the one before the "/"
      sample2 <- pair[2] #sample2 is the one after the "/"
      fst_value <- fst_results[sample1,sample2] #retrieve the Fst value for the pair of interest
      stat_df[i,fst_column_name] <- round(as.numeric(fst_value),7) #insert the Fst value in the right line and right column
    }
  }

  #[ABC ANALYSIS]

  #Create the simulated data subsets
  all_subsets <- create.subsets(stat_df)

  #Save the diversity information in files
  lapply(seq_along(all_subsets),function(i) subset.to.file(all_subsets[[i]],i)) #save the subsets in TXT files
  save(stat_df,file="stat_df.RData") #save the statistics table into a RData file

  #Create a table to retrieve the output of the ABC analysis
  abc_df <- build.abcdf(stat_df)

  #Read the observed data
  stat.obs <- read.table(paste0(folder,"/","RealDataObs.obs"),header=T) #read the observed data file
  stat.obs.list <- create.subsets.obs(stat.obs) #create subsets for the pairs of observed data

  #Perform the ABC on each subset
  selected_simul_global <- list()  #create an empty list to retrieve the best simulations
  for (i in seq_along(all_subsets)) #iterate over all the subsets (i.e. all the population pairs)
  {
    #Extract the simulated parameter values
    par.sim <- all_subsets[[i]][grep("^K1",names(all_subsets[[i]]))] #lines containing "K1"
    #Extract the simulated summary statistics
    stat.sim <- all_subsets[[i]][grep("^hapdiv",names(all_subsets[[i]]))] #lines containing "hapdiv"
    stat.sim[] <- lapply(stat.sim,as.numeric) #force the content of stat.sim to be numeric
    #Extract the observed summary statistics
    stat.obs.subset <- stat.obs.list[[i]][c(3,4)]

    #Perform the ABC analysis
    rej <- abc(target=stat.obs.subset,param=par.sim,sumstat=stat.sim,tol=0.1,method="rejection")

    #Compute the goodness of fit and the p-value
    gof_res <- gfit(target=stat.obs.subset,sumstat=stat.sim,statistic=mean,nb.replicate=nmb_simul,tol=0.1) #adjust the number of replicate to the number of simulations

    #Combine the results in a table
    retained_simulations <- all_subsets[[i]][rej$region,] #FALSE = rejected simulation ; TRUE = retained simulation
    retained_simulations$GoF <- gof_res$dist.obs  #retrieve the GoF statistic
    retained_simulations$pvalue <- summary(gof_res)$pvalue  #retrieve the p-value
    selected_simul_global[[length(selected_simul_global) + 1]] <- retained_simulations  #add the selected simulations for this subset to the list of selected simulations
  }

  #Save the information regarding the selected simulations in an RData file
  save(selected_simul_global,file="retained_simul.RData")

  #Compute the mean of the best simulations for each statistics and pair
  hapdiv_mean_1 <- sapply(selected_simul_global,function(x) mean(x[[4]])) #mean of haplotype diversity for the first sample
  hapdiv_mean_2 <- sapply(selected_simul_global,function(x) mean(x[[5]])) #mean of haplotype diversity for the second sample
  Fst_mean <- sapply(selected_simul_global,function(x) mean(x[[6]])) #mean of Fst for the pair
  GoF_mean <- sapply(selected_simul_global,function(x) mean(x[[7]]))
  pval_mean <- sapply(selected_simul_global,function(x) mean(x[[8]]))
  summary_mean <- cbind(hapdiv_mean_1,hapdiv_mean_2,Fst_mean,GoF_mean,pval_mean) #combine the means

  #Fill the ABC summary table with the values
  abc_df[2:6] <- summary_mean

  #[CONTINUITY TEST]

  test_output <- test.continuity(abc_df)

  #[CLEAN-UP]

  move.files("logs","^(slurm.*|.*\\.log)$") #move the log files in the "logs" folder
  move.files("subsets","^subset.*") #move the subset files in the "subsets" folder

  #[FINAL RETURN]

  return(test_output) #return the output of the continuity test function
}


