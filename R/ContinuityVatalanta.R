#Step 1: load selected_simul_global obtained with the function main.vatalanta

#Step 2: identify the index corresponding to the studied pairs, here: 1, 5, 8, 10.

#Save the pair data into a subset
df <- selected_simul_global[[10]]

#Extract and sort the simulated Fst values
Fst.sim.list <- selected_simul_global[[1]][[5]] #extract the values of Fst for the pair
Fst.sim.list.sorted <- sort(Fst.sim.list) #sort the Fst values

#Extract the observed Fst value
Fst.obs <- stat.obs.list[[1]][4] #retrieve the observed Fst for the pair

#Compute the p-value
pval <- round(length(Fst.sim.list.sorted[Fst.sim.list.sorted > as.numeric(Fst.obs)])/length(Fst.sim.list.sorted),3)

#Create the output data frame for the population continuity test
pop_name <- gsub("Fst","",names(df[5]))
temp_df <- data.frame(Population=pop_name,PValue=pval)
temp_df$Conclusion <- ifelse(temp_df$PValue > 0.05, "Continuity", "No continuity")
output_df <- rbind(output_df, temp_df)

# Create the visual outputs

## Fst
filename2 <- paste0(gsub("/",".",names(df[10])),".png") #create a unique file name
f <- ggplot(df, aes(x=df[[5]])) +
  geom_density(color="#1065AB")+
  labs(title = ".../...", #Update the population name here
       x = "Fst")+
  geom_vline(xintercept = as.numeric(Fst.obs),
             color = "#B31529")+
  xlim(-0.05, 0.1)+
  theme_classic()
ggsave(filename2,path = NULL) #save the plot

## POPSIZE

df$K1 <- as.numeric(as.character(df$K1))  #ensures K is numeric

# Create a unique file name
filename1 <- paste0("POPSIZE_", gsub("Fst", "", gsub("/", ".", names(df[7]))), ".png")

# Create the plot
p <- ggplot(df, aes(x = df$K1)) +
  geom_density() +
  labs(title = gsub("Fst", "", names(df[5])), x = "POPSIZE") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)  # Rotate x-axis text
  )

ggsave(filename1, plot = p,  width = 12.5, height = 9, units="cm",dpi = 300)  #save the plot
