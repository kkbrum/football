# This script is used to try the triples integer program for the football data
# This script should be run once for each stratum in 1, 3, 4, 5, 6, 7, 8
# Strata 5 and 8 are very resource intensive and will likely not finish
ist <- 1

# 0. Read in the data ---- #####################################################
install.packages("triplesmatch")
library(triplesmatch)
load("football_data_prepped_no_outcomes_r1.RData")
load("football_match_081323.RData")

# 1. Attempt integer program for triples match  ---- ###########################

# Use same number of units as in the original triples solution for that stratum
mt <- sum(m$m$nOfTreated[m$m$st == ist])
mc <- 3 * sum(m$m$st == ist) - mt

system.time({
  res <- triplesIP(z[st == ist], cost = cost[[ist]], mt = mt, mc = mc, threads = 10)
})
res$opt_info$objval
res$opt_info$objbound

# How do the objective and bound compare to the heuristic solution for that stratum?
m_st <-  triples(cost[ist], z[st == ist], st[st == ist])
m_st$obj
m_st$bound
