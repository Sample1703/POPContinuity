main.vatalanta <- function(nmb_simul,folder)
{
  #[COMPUTATION AROUND K AND ITS PRIOR]
  #Retrieve the observed value of K
  Km <- extract.Km("VAtalanta") #list of 2 containing the values of K and m
  Kuser <- Km[[2]] #value of the observed K
  #Create the prior and the list of X (nmd_simul) values of Ksim
  Kprior <- create.prior(Km) #create a list of values ranging from Kuser/10 to Kuser*10
  Klist <- generate.K.list(Kprior,nmb_simul) #create a list of X values of K taken randomly from the prior
  Kall <- c(Kuser,Klist)  #combine the values of Kuser and Ksim in a list

  #[CREATION OF AN EMPTY DATAFRAME TO RETRIEVE THE COMPUTED STATISTICS]
  stat_df <- build.stat.df("VAtalanta",nmb_simul)
  #Fill the first column of the table with the simulation ID
  stat_df[1] <- c(1:length(Kall))

  #[SIMULATION + REFORMAT ARP + STATISTICS COMPUTATION]
  #For each value of K in the Kall list
  for (i in 1:length(Kall))
  {
    #UPDATE THE VALUES OF K IN TABLE AND FILES
    #Fill the second column of the table with the value of K
    stat_df[i,2] <- Kall[i]
    #Identify the files where the value of K is present and update it with the one from the table
    replace.K.file.alt("VAtalanta",Kall[[i]])

    #LAUNCH A SIMULATION
    launch.simulation("simulation.sh",paste0(folder,".txt")) #launch the simulation using the script simulation.sh and the setting file named as the settings folder
    #Wait until the folder "GeneticsOutput is created to follow the procedure
    if (check.for.file(folder))
    {
      # Folder exists, copy its content to the working directory
      arp_file <- list.files(path=paste0(folder,"/GeneticsOutput"),pattern="\\.arp$",full.names=TRUE)
      file.copy(arp_file,".")
    }
    else
    {
      cat("GeneticsOutput folder not created.")
    }

    #CONVERT ARP FILE

    #Conversion into a dataframe
    system(paste("sbatch -p fn-mat","ARPtoFASTA.sh"))
    Sys.sleep(10)
    fasta_file <- list.files(path=".",pattern="\\.fasta$",full.names=F)
    split.fasta.file.indiv.sim(fasta_file)
    Sys.sleep(10)
    move.files("arp","^(.*\\.arp)$")

    #COMPUTE THE STATISTICS

    #Create the haplotype dataframe
    fasta_files_sim <- list.files(path=".",pattern="^[0-9].*\\.fasta$",full.names=TRUE)
    haplotype_df_sim <- build.df.vatalanta(fasta_files_sim)

    #Transform it so it can be converted into loci objects
    haplotype_loci_sim <- haplotype_df_sim
    names(haplotype_loci_sim)[names(haplotype_loci_sim) == "Population"] <- "population" #rename the population column
    haplotype_loci_sim <- haplotype_loci_sim[-1] #remove the individual ID

    haplotype_loci_sim[,-1] <- lapply(haplotype_loci_sim[,-1], function(col)
    {
      result <- gsub("([A-Za-z])([A-Za-z])", "\\1/\\2", col)
      gsub("/$","",result)
    })

    haplotype_loci_sim[,-1] <- lapply(haplotype_loci_sim[,-1], function(col)
    {
      as.factor(col)
    })

    # Create subsets for each population
    sub1910_sim <- haplotype_loci_sim[haplotype_loci_sim$population=="1910",]
    sub1950_sim <- haplotype_loci_sim[haplotype_loci_sim$population=="1950",]
    sub1960_sim <- haplotype_loci_sim[haplotype_loci_sim$population=="1960",]
    sub1970_sim <- haplotype_loci_sim[haplotype_loci_sim$population=="1970",]
    sub2022_sim <- haplotype_loci_sim[haplotype_loci_sim$population=="2022",]

    # Transform the subsets into objects of type loci
    loci1910_sim <- as.loci(sub1910_sim,ploidy=2,col.pop="population",allele.sep = "/",col.loci=2:302 )
    loci1950_sim <- as.loci(sub1950_sim,ploidy=2,col.pop="population",allele.sep = "/",col.loci=2:302 )
    loci1960_sim <- as.loci(sub1960_sim,ploidy=2,col.pop="population",allele.sep = "/",col.loci=2:302 )
    loci1970_sim <- as.loci(sub1970_sim,ploidy=2,col.pop="population",allele.sep = "/",col.loci=2:302 )
    loci2022_sim <- as.loci(sub2022_sim,ploidy=2,col.pop="population",allele.sep = "/",col.loci=2:302 )

    # Transform into genind objects
    genind1910_sim <- loci2genind(loci1910_sim,ploidy=2)
    genind1950_sim <- loci2genind(loci1950_sim,ploidy=2)
    genind1960_sim <- loci2genind(loci1960_sim,ploidy=2)
    genind1970_sim <- loci2genind(loci1970_sim,ploidy=2)
    genind2022_sim <- loci2genind(loci2022_sim,ploidy=2)

    # Compute the heterozygosity
    Hs1910_sim <- Hs(genind1910_sim)
    Hs1950_sim <- Hs(genind1950_sim)
    Hs1960_sim <- Hs(genind1960_sim)
    Hs1970_sim <- Hs(genind1970_sim)
    Hs2022_sim <- Hs(genind2022_sim)
    Hs_results <- list(Hs1910_sim[[1]],Hs1950_sim[[1]],Hs1960_sim[[1]],Hs1970_sim[[1]],Hs2022_sim[[1]])

    # Convert the full simulated dataframe into a loci object
    haplotype_loci_sim <- as.loci(haplotype_loci_sim,ploidy=2,col.pop="population",allele.sep = "/",col.loci=2:302 )
    # Convert the loci object into a genind object
    genind_sim <- loci2genind(haplotype_loci_sim,ploidy=2)
    # Compute the Fsts
    fst_results_va <- genet.dist(genind_sim,method="Nei87")

    #Retrieve information on the samples
    num_samples <- length(unique(haplotype_df_sim$Population))
    sample_names <- unique(haplotype_df_sim$Population)

    #Add colnames to stat_df
    sample_combinations <- combn(sample_names,2,function(x) paste0("Fst",paste(x,collapse="/")),simplify=T) #create a combination of : "Fst",sample name,"/",sample name
    #print(sample_combinations)
    colnames(stat_df) <- c("SimID","K1","MutationRate",paste0("intra",sample_names),sample_combinations) #rename the columns

    #Fill the table with the intrapop information
    stat_df[i,4:(as.numeric(num_samples)+3)] <- unlist(round(as.numeric(Hs_results),7)) #fill the columns 3 > num_samples+1 with the intrapop values

    #Fill the table with the Fst information
    stat_df[i,9:(8 + as.numeric(num_samples * (num_samples - 1) / 2))] <- round(as.numeric(fst_results_va),7)

    move.files("fasta","^(.*\\.fasta)$")
  }

  stat_df <- stat_df[, colSums(!is.na(stat_df)) > 0] #remove the empty columns

  #[ABC ANALYSIS]
  #Create the data subsets
  all_subsets <- create.subsets.alt(stat_df)
  #Save the subsets in TXT files
  lapply(seq_along(all_subsets),function(i) subset.to.file(all_subsets[[i]],i))
  #Save the statistics dataframe into a RData file
  save(stat_df,file="stat_df_K5000_9.RData")
  #Create the summary table
  abc_df <- build.abcdf(stat_df)

  #For each subset we want to select the best simulations

  #Read the observed data
  stat.obs <- read.table(paste0(folder,"/","RealDataObs.obs"),header=T) #read the observed data file
  stat.obs.list <- create.subsets.obs(stat.obs) #create subsets for the pairs of observed data

  #Perform the ABC on each subset
  selected_simul_global <- list()  # Create an empty list to retrieve the best simulations
  for (i in seq_along(all_subsets)) # Use seq_along() to iterate over all_subsets
  {
    #Extract the simulated parameter values
    par.sim <- all_subsets[[i]][grep("^K1",names(all_subsets[[i]]))] #lines containing "K1"
    #Extract the simulated summary statistics
    stat.sim <- all_subsets[[i]][grep("^intra",names(all_subsets[[i]]))] #lines containing "intra"
    stat.sim[] <- lapply(stat.sim,as.numeric) #force the content of stat.sim as numerical values
    stat.obs.subset <- stat.obs.list[[i]][c(2,3)]

    #Perform the ABC analysis
    rej <- abc(target=stat.obs.subset[1],param=par.sim,sumstat=stat.sim[1],tol=0.1,method="rejection")
    #Retrieve the best simulations
    ##Compute the goodness of fit and the p-value
    gof_res <- gfit(target=stat.obs.subset[1],sumstat=stat.sim[1],statistic=mean,nb.replicate=10,tol=0.1) #adjust the number of replicate to the number of simulations
    ##Combine the results in a table
    retained_simulations <- all_subsets[[i]][rej$region,] # FALSE = rejected simulation ; TRUE = retained simulation
    retained_simulations$GoF <- gof_res$dist.obs  # Retrieve the GoF statistic
    retained_simulations$pvalue <- summary(gof_res)$pvalue  # Retrieve the p-value
    selected_simul_global[[length(selected_simul_global) + 1]] <- retained_simulations  # Add the selected simulations for this subset to the list of selected simulations
  }

  #Save the informations of the selected simulations
  save(selected_simul_global,file="retained_simul_K5000_9.RData")

  #Compute the mean of the best simulations for each statistics and pair
  hapdiv_mean_1 <- sapply(selected_simul_global,function(x) mean(x[[3]])) #mean of intra for the first sample
  hapdiv_mean_2 <- sapply(selected_simul_global,function(x) mean(x[[4]])) #mean of intra for the second sample
  Fst_mean <- sapply(selected_simul_global,function(x) mean(x[[5]])) #mean of Fst for the pair
  GoF_mean <- sapply(selected_simul_global,function(x) mean(x[[6]]))
  pval_mean <- sapply(selected_simul_global,function(x) mean(x[[7]]))
  summary_mean <- cbind(hapdiv_mean_1,hapdiv_mean_2,Fst_mean,GoF_mean,pval_mean) #combine the means

  #Fill the abc summary table with the values
  abc_df[2:6] <- summary_mean

  return(abc_df)
}
