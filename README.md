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

##### Computations around K and its prior

- Retrieve the value of K from the user dynamic_K file (only 1 vegetation file, no change, original value 5,000 to have a range between 500 and 50,000)





