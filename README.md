# POPContinuity
Testing for continuity in structured populations with a generic R tool.

This project aims to develop an R package, called POPContinuity, that performs a structured population continuity test, investigating whether there has been a detectable population replacement in the time and area under study. The specificity of this approach is that it takes into account the continuous gene flow between the population of interest and neighboring populations, a factor that has proved important when investigating population continuity over time. Its main features are, first, the simulation of genetic data of different types (SNP, microsatellites, or sequences) under a null hypothesis of population continuity, and second, the comparison with observed data while accounting for both spatial and temporal dimensions of the population. POPcontinuity relies on the Approximate Bayesian Computation (ABC) analysis to optimize the goodness-of-fit of the simulations, based on intra-population genetic diversity, while the continuity test is based on inter-population genetic divergence. Integrating simulations and analysis into a single, easy-to-use R package simplifies the user approach to population continuity testing. POPcontinuity also provides various complete summary files for detailed results, making it a flexible tool, adaptable to various research contexts when molecular data in temporal series is available. 
To illustrate the usefulness of our package and to validate its results, we applied it to two case studies from different research fields: the evolutionary history of modern human populations in Europe and the conservation of the butterfly Vanessa atalanta in Switzerland. 

## Specific methodology for the application on human evolution

The data files concerning the application on human evolution are in the folder "inst/extdata/CentralEurope"
Below is a description of the files located in this folder. 

**SPLATCHE3_CentralEuropeSimulation**: Folder containing all the settings file for SPLATCHE3

**SPLATCHE3_CentralEuropeSimulation.txt**: Settings file with the information required to launch the simulation

**simulations.sh**: Bash script to launch the simulations

**SPLATCHE4-Linux-64b**: SPLATCHE3 executable

### How to use this application

The function main.haploid() located in "R/MainHaploid" requires two inputs : 
- The number of simulations to perform
- The name of the folder containing the settings

Use the function as follow : **main.haploid(10000,"SPLATCHE3_CentralEuropeSimulation")**

