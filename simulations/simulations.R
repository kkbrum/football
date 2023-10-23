# This script is used for simulations of the triples design for the football appendix


# 0. Preparation ---- ################################

## Packages ----
install.packages("triplesmatch")
library(triplesmatch)
library(pracma)
library(MASS)
library(sensitivityfull)
library(rrelaxiv)

## Settings ----
set.seed(128, kind = "Mersenne-Twister")
ncov <- 6
n_t <- 200
pi_t_list <- c(1/3, 2/5, 1/2) # Proportion of observations that are treated
delta_list <- c(0.5, 1) # Bias in covariates
n_settings <- length(pi_t_list) * length(delta_list)
trim <- 3
gamma_list <- c(1, 2.5, 3, 3.25, 3.5, 3.75, 4, 4.25)
nsim <- 200 # Run each setting 200 times
n_type <- 5 # How many different types of matches are we investigating?

## Dataframes to hold results ----

### For power ----
power_res <- data.frame(matrix(0, nrow = n_settings * n_type, ncol = length(gamma_list) + 6,
                               dimnames = list(NULL, c("bias", "nt", "nc", "match_type",
                                                       "n_t_used", "n_c_used",
                                                       paste0("Gamma=", gamma_list)))))
power_res$match_type <- rep(c("Pair", "Pair (no crossing)", "Pair (reduced)", "Triple",
                              "Full (no crossing)"), n_settings)
ncol_p <- ncol(power_res)

### For match quality ----
res <- data.frame( matrix(0, nrow = n_settings * n_type, ncol = 20,
                          dimnames = list(NULL,
                                          c("bias", "nt", "nc",
                                            "before_cost_only_within", "before_cost_all",
                                            "match_type",
                                            "n_crosses", "n_t_used", "n_c_used", "n12", "n11", "n21",
                                            "cost_subset", "cost_total_unpen",
                                            "cost_avg_unpen",  "cost_avg_pen",
                                            "cost_med", "cost_q75", "cost_q90", "cost_q95"
                                          ))))
res$match_type <- rep(c("Pair", "Pair", "Pair", "Triple", "Full"), n_settings)


# 1. Simulations ---- ##########################################################

for (sim in 1:nsim) {
  print(paste0("Simulation number ", sim))
  res_row <- 1
  for(delta in delta_list) {
    print(paste0("Covariate bias of ", delta))
    for(pi_t in pi_t_list) {
      print(paste0("Treated proportion of ", pi_t))

      res[res_row + 0:(n_type - 1), "bias"] <-  delta
      power_res[res_row + 0:(n_type - 1), "bias"] <- delta

      ## Data generation ----

      nobs <- round(n_t / pi_t) # Number of total individuals
      n_c <- nobs - n_t # Number of control individuals
      res[res_row + 0:(n_type - 1), "nt"] <- n_t
      res[res_row + 0:(n_type - 1), "nc"] <- n_c
      power_res[res_row + 0:(n_type - 1), "nt"] <- n_t
      power_res[res_row + 0:(n_type - 1), "nc"] <- n_c
      z <- c(rep(1, n_t), rep(0, n_c)) # Treated indicator

      ### Generate covariate data ----

      # Spread the covariate bias randomly over all covariates
      mu_t <- as.vector(c(delta, rep(0, ncov-1)) %*% randortho(ncov, type = "orthonormal"))
      mu_c <- rep(0, ncov)
      id <- diag(ncov)
      X <- rbind(mvrnorm(n_t, mu_t, id), mvrnorm(n_c, mu_c, id))

      # Estimate propensity score from data and add as covariate
      ps <- glm(z ~ X, family = binomial)$fitted.values
      # X <- cbind(X,ps) # let's not add the PS to X for the moment so it is excluded from the distance matrix
      # 5 propensity score strata to match within
      ps_st <- cut(ps, c(0, quantile(ps, 1/5 * 1:4), 1), labels = 1:5)

      ### Create cost matrices ----

      # Use the covs and ps for the distance matrix
      # All distances including across strata, with no penalty for crossing strata
      cost_unpen <-  dist_mahal(X = X, st = factor(rep(1, length(z))), z = z, rank_based = FALSE)
      # Distances within strata
      cost_st <- dist_mahal(X = X, st = ps_st, z = z, rank_based = FALSE)

      cost <- cost_unpen
      # Create a penalized cost with a penalty of 1000 for crossing strata
      cost[[1]] <- cost[[1]] + outer(ps_st[z==1], ps_st[z==0], "!=") * 1000
      # Create a penalized cost with infinite penalty for crossing
      cost_inf <- cost_unpen
      cost_inf[[1]][outer(ps_st[z==1], ps_st[z==0], "!=")] <- Inf

      ### Generate outcome ----

      # N(0, 1/2) for controls, N(1/2, 1/2) for treated. Thus any t-c difference is N(1/2, 1)
      # Additive treatment effect tau,
      #     but no bias from covariates so that quality of match doesn't skew results
      tau <- 1/2
      Y <- c(rnorm(n_t, tau, sqrt(1/2)), rnorm(n_c, 0, sqrt(1/2)))
      names(Y) <- 1:length(Y)

      ### Calculate avg cost before matching ----

      # Just within strata
      before_cost <- 0
      before_n_cost <- 0
      for (ist in as.numeric(levels(ps_st))) {
        before_cost <- before_cost + sum(cost_st[[ist]])
        before_n_cost <- before_n_cost + length(cost_st[[ist]])
      }
      res[res_row + 0:(n_type - 1), "before_cost_only_within"] <- before_cost / before_n_cost

      # Across strata too, but not penalizing for stratum crossings
      res[res_row + 0:(n_type - 1), "before_cost_all"] <- mean(cost_unpen[[1]])

      ## Matching ----

      ### Run triples match first ----
      triplem <-  suppressWarnings(triples(cost = cost_st, z = z, st = ps_st)$m)
      triplem_info <- triplesmatch:::process_triples_match(m = triplem, gamma_list = gamma_list, Y = Y, trim = trim)

      res[res_row + 3, 7:20] <- res[res_row + 3, 7:20] + triplem_info$quality

      power_res[res_row + 3, 5:6] <- power_res[res_row + 3, 5:6] +
        triplem_info$quality[c("n_t_used", "n_c_used")]
      power_res[res_row + 3, 7:ncol_p] <- power_res[res_row + 3, 7:ncol_p] + triplem_info$reject

      # How many treated units to drop in the pair and full match if desired to have equal treated counts
      treated_to_drop <- 200 - triplem_info$quality['n_t_used']

      ### Run pair match ----

      #### Pair match using only the number of treated units matched in the triples match ----

      if (treated_to_drop > 0) {
        # Add extra control units: any treated units matched to these controls will be dropped
        cost_w_extra_controls <- cbind(cost[[1]],
                                       matrix(0, nrow = nrow(cost[[1]]),
                                              ncol = treated_to_drop,
                                              dimnames =
                                                list(NULL, paste0("drop", 1:treated_to_drop))))
      } else {
        cost_w_extra_controls <- cost[[1]]
      }
      suppressWarnings(pairm <-  optmatch::pairmatch(x = cost_w_extra_controls, controls = 1))

      # Drop any units matched to fake controls
      pairm[pairm %in% pairm[paste0('drop', 1:treated_to_drop)]] <- NA
      pairm <- pairm[1:(n_t+n_c)]
      pairm <- droplevels(pairm)

      pairm_info <- triplesmatch:::process_optmatch(m = pairm, cost_unpen = cost_unpen[[1]],
                                     cost_pen = cost_inf[[1]], gamma_list = gamma_list,
                                     Y = Y, trim = trim)

      res[res_row + 2, 7:20] <- res[res_row + 2,  7:20] + pairm_info$quality

      power_res[res_row + 2, 5:6] <- power_res[res_row + 2, 5:6] +
        pairm_info$quality[c("n_t_used", "n_c_used")]
      power_res[res_row + 2, 7:ncol_p] <- power_res[res_row + 2, 7:ncol_p] + pairm_info$reject


      #### Pair match using all treated units, allowing stratum crossings ----

      if (triplem_info$quality['n_t_used'] != n_t) {
        # If we did not already intend to use all treated units, run a new match
        suppressWarnings(pairm_all_units <-  optmatch::pairmatch(x = cost[[1]], controls = 1))
        pairm_all_units_info <-
          triplesmatch:::process_optmatch(m = pairm_all_units, cost_unpen = cost_unpen[[1]],
                           cost_pen = cost_inf[[1]], gamma_list = gamma_list,
                           Y = Y, trim = trim)

      } else {
        # If the previous match was already supposed to use all treated units
        # due to triples not dropping any, just use the same match
        pairm_all_units_info <- pairm_info
      }

      res[res_row + 0, 7:20] <- res[res_row + 0, 7:20] +
        pairm_all_units_info$quality

      power_res[res_row + 0, 5:6] <- power_res[res_row + 0, 5:6] +
        pairm_all_units_info$quality[c("n_t_used", "n_c_used")]
      power_res[res_row + 0, 7:ncol_p] <- power_res[res_row + 0, 7:ncol_p] +
        pairm_all_units_info$reject


      #### Pair match, no stratum crossings ----

      if (pairm_info$quality['n_t_used'] == n_t &
          pairm_info$quality['n_crosses'] == 0) {
        # If the first match used all treated units and did not cross strata, just use that
        pairm_no_cross_info <- pairm_info
      } else if (pairm_all_units_info$quality['n_crosses'] == 0) {
        # If the second match using all treated units did not cross strata, use that
        pairm_no_cross_info <- pairm_all_units_info
      } else {
        # Create new match if needed
        # Run separately for each stratum and then combine across strata
        pairm_no_cross_st <- list()
        for (ist in levels(ps_st)) {
          suppressWarnings(
            pairm_no_cross_st[[ist]] <- optmatch::pairmatch(x = cost_st[[as.numeric(ist)]], controls = 1))
        }
        pairm_no_cross <- c(pairm_no_cross_st[[1]], pairm_no_cross_st[[2]],
                            pairm_no_cross_st[[3]], pairm_no_cross_st[[4]],
                            pairm_no_cross_st[[5]])

        pairm_no_cross_info <- triplesmatch:::process_optmatch(m = pairm_no_cross,
                                                cost_unpen = cost_unpen[[1]],
                                                cost_pen = cost_inf[[1]],
                                                gamma_list = gamma_list,
                                                Y = Y, trim = trim)
      }

      res[res_row + 1, 7:20] <-  res[res_row + 1, 7:20] +
        pairm_no_cross_info$quality

      power_res[res_row + 1, 5:6] <- power_res[res_row + 1, 5:6] +
        pairm_no_cross_info$quality[c("n_t_used", "n_c_used")]
      power_res[res_row + 1, 7:ncol_p] <- power_res[res_row + 1, 7:ncol_p] +
        pairm_no_cross_info$reject


      ### Run full match using all treated units w/ no stratum crossings ----

      tab <- table(z, ps_st)
      c_drop <- (pmax(tab[1, ] - 2 * tab[2, ], 0))
      t_drop <- (pmax(tab[2, ] - 2 * tab[1, ], 0))

      fullm_no_cross_st <- list()
      for (ist in levels(ps_st)) {
        minc <- 1/2
        maxc <- 2
        omit_frac <- 0
        if (c_drop[ist] > 0) {
          minc <- 2
          omit_frac <- c_drop[ist]/tab[1, ist]
        } else if (t_drop[ist] > 0) {
          maxc <- 1/2
          omit_frac <- -t_drop[ist]/tab[2, ist]
        }

        suppressWarnings(
          fullm_no_cross_st[[ist]] <- optmatch::fullmatch(x = cost_st[[as.numeric(ist)]],
                                                          min.controls = minc, max.controls = maxc,
                                                          omit.fraction = omit_frac))
      }
      fullm_no_cross <- c(fullm_no_cross_st[[1]], fullm_no_cross_st[[2]],
                          fullm_no_cross_st[[3]], fullm_no_cross_st[[4]],
                          fullm_no_cross_st[[5]])

      fullm_no_cross_info <- triplesmatch:::process_optmatch(m = fullm_no_cross,
                                              cost_unpen = cost_unpen[[1]],
                                              cost_pen = cost_inf[[1]],
                                              gamma_list = gamma_list,
                                              Y = Y, trim = trim)

      if (fullm_no_cross_info$quality['n_c_used'] != n_c - sum(c_drop)) {
        stop("Wrong number of controls dropped")}
      if (fullm_no_cross_info$quality['n_t_used'] != n_t - sum(t_drop)) {
        stop("Wrong number of treated units dropped")}

      res[res_row + 4, 7:20] <- res[res_row + 4, 7:20]  +
        fullm_no_cross_info$quality

      power_res[res_row + 4, 5:6] <- power_res[res_row + 4, 5:6] +
        fullm_no_cross_info$quality[c("n_t_used", "n_c_used")]
      power_res[res_row + 4, 7:ncol_p] <- power_res[res_row + 4, 7:ncol_p] +
        fullm_no_cross_info$reject

      res_row <- res_row + 5

    }
  }
}

# Average results across the simulations
res[, 7:20] <- res[, 7:20] / nsim
power_res[, 5:ncol_p] <- power_res[, 5:ncol_p] / nsim

save.image(file = "simulations.RData")
# load("simulations.RData")

# 2. Tables ---- ###############################################################
library(xtable)

## Match quality table ----

res[(1:nrow(res)) %% n_type != 1, 2:5] <- NA
suppressWarnings(res$n_crosses[!is.na(as.numeric(res$n_crosses))] <-
                   round(as.numeric(res$n_crosses[!is.na(as.numeric(res$n_crosses))]), 2))
res[(1:nrow(res)) %% n_type %in% c(2,4,0), "n_crosses"] <- "\\times"

res_reduced <- res
res_reduced[(1:nrow(res_reduced)) %% 5 %in% 1:3, "match_type"] <- rep(paste0("Pair ", 1:3), n_settings)

res_reduced[, 8:9] <- round(res_reduced[, 8:9], 2)
for (i in c(8:9)) {
  for (j in 1:12) {
    ind <- (j-1)*5 + 1:5
    temp <- res_reduced[ind, i]
    res_reduced[ind[which(temp == max(temp))], i] <-
      paste0("\\textbf\\{", res_reduced[ind[which(temp == max(temp))], i], "\\}")
  }
}

res_reduced[, 14:20] <- round(res_reduced[, 14:20], 2)
for (i in c(14:20)) {
  for (j in 1:12) {
    ind <- (j-1)*5 + 1:5
    temp <- res_reduced[ind, i]
    res_reduced[ind[which(temp == min(temp))], i] <-
      paste0("\\textbf\\{", res_reduced[ind[which(temp == min(temp))], i], "\\}")
  }
}

for (delta in c(0.5, 1)) {
  print.xtable(xtable(res_reduced[res_reduced$bias == delta, c(-1, -13, -14)],
                      caption = paste0("Simulated match quality results for bias of ", delta),
                      digits = c(0, 0, 0, rep(2, 15))),
               santize.text.function = identity,
               include.rownames = FALSE, hline.after = c(-1, 0, 5, 10, 15))
}

## Simulated power table ----

power_res[(1:nrow(res)) %% n_type != 1, 2:3] <- NA
power_res_reduced <- power_res
power_res_reduced <- power_res_reduced[(1:nrow(power_res_reduced)) %%
                                         (n_type * length(pi_t_list)) %in% 6:10, ]
power_res_reduced[(1:nrow(power_res_reduced)) %% 5 %in% 1:3, "match_type"] <-
  rep(paste0("Pair ", 1:3), length(delta_list))
power_res_reduced[(1:nrow(power_res_reduced)) %% 5 %in% 0, "match_type"] <- "Full"

power_res_reduced[, 7:14] <- round(power_res_reduced[, 7:14], 2)
for (i in c(7:14)) {
  for (j in 1:4) {
    ind <- (j-1)*5 + 1:5
    temp <- power_res_reduced[ind, i]
    power_res_reduced[ind[which(temp == max(temp))], i] <-
      paste0("\\textbf\\{", power_res_reduced[ind[which(temp == max(temp))], i], "\\}")
  }
}

for (delta in c(0.5, 1)) {
  print(xtable(power_res_reduced[power_res_reduced$bias == delta, -1],
               caption = paste0("Simulated power results for bias of ", delta),
               digits = c(0, 0, 0, rep(2, 11))),
        santize.text.function = identity,
        include.rownames = FALSE, hline.after = c(-1, 0, 5))
}
