# football
Reproducibility materials for "A New Design for Observational Studies Applied to the Study of the Effects of High School Football on Cognition Late in Life."

# Data availability
Public data are available at https://researchers.wls.wisc.edu/data/survey-data/.
Our analyses require access to the yearbook data, which is private.
The private version of the data must be requested in order to fully reproduce our analyses; instructions for obtaining the private data are at the same link.

# R package
Our analyses require the 'triplesmatch' R package, developed for this project. This package can be installed and loaded by running the following in R:

```r
install.packages("triplesmatch")
library(triplesmatch)
```

# Organization

The code folder contains four files. They should be run in the following order:

1. football_prep.R: This file loads in the data obtained from the Wisconsin Longitudinal Study and implements the eligibility and exclusion criteria for our study.
2. football_triples_matching.R: This file loads in the prepared data and implements the triples match. Additionally generates figures 1-2 and tables 1-3 of the paper. Finally, creates coherent ranks for matched individuals.
3. football_analysis.R: This file creates figure 4 and implements the final analysis of the matched individuals. 
4. football_triples_IP.R: This file attempts to conduct the triples match using an integer program; this is not relevant to the main analysis and is discussed more in the appendix of the paper. Recreates table 3 of the appendix.

The simulations folder contains two files.

1. simulations.R: This file creates the simulations that are discussed in the appendix of the paper. Recreates tables 1-2 of the appendix.
2. simulations.RData: This file contains the data generated by the simulations (since they take a while to run, this is useful if recreating the tables in the appendix).
