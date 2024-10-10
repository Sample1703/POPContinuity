#' Test continuity of populations based on Fst values
#'
#' This function assesses the continuity of populations by comparing simulated Fst values with observed Fst values. It generates a dataframe with the results, including the population pair, the p-value of the test, and a conclusion indicating "Continuity" or "No continuity". Additionally, it saves visual representations of population sizes and Fst values.
#'
#' @param selected_simul_global A list containing dataframes of simulated statistics for different population pairs. Each dataframe should have the Fst values stored in the 7th column.
#' @param stat.obs.list A list of observed statistics, with each entry corresponding to a population pair. The observed Fst value should be the fourth element in each list entry.
#'
#' @return A dataframe with columns `Population`, `PValue`, and `Conclusion`. The `Population` column indicates the population pair, the `PValue` column provides the p-value from the comparison, and the `Conclusion` column states whether continuity is observed (`"Continuity"`) or not (`"No continuity"`).
#'
#' @details
#' The function calculates p-values by comparing the observed Fst values to the distribution of simulated Fst values. A p-value greater than 0.05 suggests continuity between the populations, while a p-value of 0.05 or less suggests no continuity.
#' The function also generates and saves density plots of population sizes and Fst values for each pair. The plots are saved in the current working directory, and are moved into folders named `"POPSIZE"` and `"FST"` respectively.
#'
#' @examples
#' results <- test.continuity(selected_simul_global, stat.obs.list)
#'
#' @export
test.continuity <- function(selected_simul_global,stat.obs.list)
{
  output_df <- data.frame(Population=character(),Pvalue=numeric(),Conclusion=character())

  for (i in 1:length(selected_simul_global))
  {
    df <- selected_simul_global[[i]]

    #Simulated Fst
    Fst.sim.list <- selected_simul_global[[i]][[7]] #extract the values of Fst for the pair
    Fst.sim.list.sorted <- sort(Fst.sim.list) #sort the Fst values

    #Observed Fst
    Fst.obs <- stat.obs.list[[i]][4] #retrieve the observed Fst for the pair

    #Comparison
    pval <- round(length(Fst.sim.list.sorted[Fst.sim.list.sorted > as.numeric(Fst.obs)])/length(Fst.sim.list.sorted),3)
    pop_name <- gsub("Fst","",names(df[7]))
    temp_df <- data.frame(Population=pop_name,PValue=pval)
    temp_df$Conclusion <- ifelse(temp_df$PValue > 0.05, "Continuity", "No continuity")
    output_df <- rbind(output_df, temp_df)

    # Visual representations

    ## POPSIZE
    filename1 <- paste0("POPSIZE_",gsub("Fst","",gsub("/",".",names(df[7]))),".png") #create a unique file name
    p <- ggplot(df, aes(x=df$K1))+
      geom_density()+
      labs(title=gsub("Fst","",names(df[7])),x="POPSIZE")
    ggsave(filename1,path = NULL)

    ## Fst
    filename2 <- paste0(gsub("/",".",names(df[7])),".png") #create a unique file name
    f <- ggplot(df, aes(x=df[[7]])) +
      geom_density()+
      labs(title=gsub("Fst","",names(df[7])),x="Fst")+
      geom_vline(xintercept = as.numeric(Fst.obs),color="red")
    ggsave(filename2,path = NULL)
  }

  move.files("POPSIZE","^POPSIZE.*") #move the popsize plots in the "POPSIZE" folder
  move.files("FST","^Fst*") #move the fst plots in the "FST" folder
  return(output_df)
}

test.continuity.alt <- function(selected_simul_global)
{
  output_df <- data.frame(Population=character(),Pvalue=numeric(),Conclusion=character())

  for (i in 1:length(selected_simul_global))
  {
    df <- selected_simul_global[[i]]

    #Simulated Fst
    Fst.sim.list <- df[[5]] #extract the values of Fst for the pair
    Fst.sim.list.sorted <- sort(Fst.sim.list) #sort the Fst values

    #Observed Fst
    Fst.obs <- stat.obs.list[[i]][4] #retrieve the observed Fst for the pair

    #Comparison
    pval <- round(length(Fst.sim.list.sorted[Fst.sim.list.sorted > as.numeric(Fst.obs)])/length(Fst.sim.list.sorted),3)
    pop_name <- gsub("Fst","",names(df[5]))
    temp_df <- data.frame(Population=pop_name,PValue=pval)
    temp_df$Conclusion <- ifelse(temp_df$PValue > 0.05, "Continuity", "No continuity")

    output_df <- rbind(output_df,temp_df)
  }
  return(output_df)
}
