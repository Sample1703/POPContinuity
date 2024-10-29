# POPContinuity
Testing for continuity in structured populations with a generic R tool.

<p align="justify">
This project aims to develop an R package, called POPContinuity, that performs a structured population continuity test, investigating whether there has been a detectable population replacement in the time and area under study. The specificity of this approach is that it takes into account the continuous gene flow between the population of interest and neighboring populations, a factor that has proved important when investigating population continuity over time. Its main features are, first, the simulation of genetic data of different types (SNP, microsatellites, or sequences) under a null hypothesis of population continuity, and second, the comparison with observed data while accounting for both spatial and temporal dimensions of the population. POPcontinuity relies on the Approximate Bayesian Computation (ABC) analysis to optimize the goodness-of-fit of the simulations, based on intra-population genetic diversity, while the continuity test is based on inter-population genetic divergence. Integrating simulations and analysis into a single, easy-to-use R package simplifies the user approach to population continuity testing. POPcontinuity also provides various complete summary files for detailed results, making it a flexible tool, adaptable to various research contexts when molecular data in temporal series is available. 
</p>

## Specific methodology for the application on Vanessa Atalanta

The data files concerning the application on Vanessa Atalanta are in the folder "inst/extdata/VAtalanta"
Below is a description of the files located in this folder. 

**VAtalanta**: Folder containing all the settings file for SPLATCHE3

**VAtalanta.txt**: Settings file with the information required to launch the simulation

**simulations.sh**: Bash script to launch the simulations

**ARPtoFASTA.sh**: Bash script to convert an ARP file into a FASTA file

**SPLATCHE4-Linux-64b**: SPLATCHE3 executable

### How to perform the analysis

#### Main function 

**Computations around K and its prior**

- Retrieve the value of K from the user dynamic_K file (number of vegetation file: 1, change: no, original value: 5,000)
- Build a prior of K values ranging from K/10 to K*10 (from 500 to 50,000)
- Create a list of random K values of length equal to the number of simulations to perform (g.e. if 10,000 are performed, extract 10,000 random values of K)
- Build an empty data frame designed to store the diversity statistics values
- In the first column of the data frame, insert the simulation ID

**In a for loop**

##### Generate and convert the .ARP file 

- In the second column of the data frame, insert the random value of K (new K value)
- Replace the old K value in the vegetation with the new K value
- Launch 1 simulation
- If the folder "GeneticsOutput" is created, copy the .ARP inside it into the working directory
- Convert the .ARP file into a FASTA file
- Move the .ARP file in a folder called "arp"

##### Compute the diversity statistics

- Step 1: Create a haplotype data frame containing the simulated genetic data
- Step 2: Perform some small reformatting of the data frame so it can be converted into a loci object

###### To compute the expected heterozygosity

- For each sample, create a dedicated subset (g.e. one for 1910, one for 1950, and so on)
- Convert each subset into a loci object
- Convert each loci object into a genind object
- For each genind object, compute the expected heterozygosity

###### To compute the Fst

- Using the output of step 2, convert the data frame into a loci object
- Convert the loci object into a genind object
- Compute the Fst

Once computed, the diversity statistics are stored in the statistics data frame created earlier. 

**ABC analysis**

- Create a subset of the data frame containing the diversity statistics for each sample pair
- Save the subsets into .TXT files
- Save the statistics data frame into an RData file
- Create a new empty data frame to store the output of the ABC analysis
- Read the observed statistics, located in a file called "RealDataObs.obs" and create a subset of each sample pair 

##### For each subset (for loop)

- Step 1: Perform the ABC analysis using the intra-population diversity values of the two samples in the pair
- Step 2: Save the statistics of the simulations retained by the ABC analysis into an RData file
- Step 3: Save the results of the ABC analysis into the previously created empty data frame

The main function returns the output of the ABC analysis, a data frame containing :

- The name of the studied pair
- The intra-population diversity value for the first sample of the pair
- The intra-population diversity value for the second sample of the pair
- The Fst value of the sample pair
- The goodness-of-fit
- Its associated p-value

This process generates many files (log files, subsets, FASTA files, etc). Some code lines are available to move them into dedicated folders to free up the working directory. 

```
abc_df_K5000 <- main.vatalanta(10000,"VAtalanta")
```

The function outputs the ABC results when the starting K value is 5,000; 10,000 simulations are performed and the settings folder to use is "VAtalanta".

#### Continuity test 

Using one of the outputs of the ABC analysis, it is possible to identify the pairs for which the model is able to reproduce the observed intra-population diversity. 

During the ABC analysis, the statistics of the simulations retained by the ABC analysis are saved into an RData file (cf. Step 2 above). 

- Import the RDatafile containing the statistics of the retained simulations
- Retrieve the selected simulations for the first sample pair
- Extract the simulated Fst values
- Extract the observed Fst value
- Compute the p-value for the sample pair using the previously extracted Fst values
- Build the output data frame containing: the name of the sample pair, the associated p-value, the conclusion of the test
- Create and save the visual output for the continuity test (simulated Fst distribution and observed Fst value)
- Create and save the visual output for the carrying capacity value

Cf. ```ContinuityVatalanta.R```
The content of this file should be put into a function and inserted into the main function ```main.vatalanta()```.


