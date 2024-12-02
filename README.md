# A Two-Stage Bayesian Small Area Estimation Approach for Proportions

This repository contains `R` and `stan` code for the analyses described in the manuscript "A Two-Stage Bayesian Small Area Estimation Approach for Proportions" by James Hogg, Jessica Cameron, Susanna Cramb, Peter Baade, and Kerrie Mengersen. The published article is available [here](https://onlinelibrary.wiley.com/doi/full/10.1111/insr.12572). 

## Simulation Study

In the `simulation_study` folder is the code to reproduce the simulation study. Note that the code was run using a HPC system at the Queensland University of Technology. The `figures_tables.R` is used to reproduce the tables and figures in the paper from the results found [here](https://drive.google.com/file/d/1_KBia2SH6IpqBiHLqDffwUaIKQieBDx5/view?usp=sharing). Please contact the corresponding author for details on how to use the simulation study code. 

## Case Study

In the `case_study` folder is the code to reproduce the case study plots using the model results found [here](https://drive.google.com/file/d/1_KBia2SH6IpqBiHLqDffwUaIKQieBDx5/view?usp=sharing). The unit record National Health Survey data cannot be released. The `summary_fit.R` file is used to summarize the raw posterior draws from the model fits. `figures_tables.R` is used to create the plots and tables. Note that `functions_ALL.R` includes many user-written functions required to reproduce the case study plots. 



