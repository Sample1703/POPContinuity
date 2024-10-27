# POPContinuity
Testing for continuity in structured populations with a generic R tool.

This project aims to develop an R package, called POPContinuity, that performs a structured population continuity test, investigating whether there has been a detectable population replacement in the time and area under study. The specificity of this approach is that it takes into account the continuous gene flow between the population of interest and neighboring populations, a factor that has proved important when investigating population continuity over time. Its main features are, first, the simulation of genetic data of different types (SNP, microsatellites, or sequences) under a null hypothesis of population continuity, and second, the comparison with observed data while accounting for both spatial and temporal dimensions of the population. POPcontinuity relies on the Approximate Bayesian Computation (ABC) analysis to optimize the goodness-of-fit of the simulations, based on intra-population genetic diversity, while the continuity test is based on inter-population genetic divergence. Integrating simulations and analysis into a single, easy-to-use R package simplifies the user approach to population continuity testing. POPcontinuity also provides various complete summary files for detailed results, making it a flexible tool, adaptable to various research contexts when molecular data in temporal series is available. 
To illustrate the usefulness of our package and to validate its results, we applied it to two case studies from different research fields: the evolutionary history of modern human populations in Europe and the conservation of the butterfly Vanessa atalanta in Switzerland. 

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

#### Generate and convert the .ARP file 

- In the second column of the data frame, insert the random value of K (new K value)
- Replace the old K value in the vegetation with the new K value
- Launch 1 simulation
- If the folder "GeneticsOutput" is created, copy the .ARP inside it into the working directory
- Convert the .ARP file into a FASTA file
- Move the .ARP file in a folder called "arp"

#### Compute the diversity statistics

- Step 1: Create a haplotype data frame containing the simulated genetic data
- Step 2: Perform some small reformatting of the data frame so it can be converted into a loci object

##### To compute the expected heterozygosity

- For each sample, create a dedicated subset (g.e. one for 1910, one for 1950, and so on)
- Convert each subset into a loci object
- Convert each loci object into a genind object
- For each genind object, compute the expected heterozygosity

##### To compute the Fst

- Using the output of step 2, convert the data frame into a loci object
- Convert the loci object into a genind object
- Compute the Fst

One computed diversity statistics are stored in the statistics data frame created earlier. 




