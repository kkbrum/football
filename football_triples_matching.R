# This script is used to create the triples match (and create relevant figures and tables)
# Also creates aberrant and coherent ranks
#      (and creates relevant figures and tables of their marginal distributions)
# This is script 2 of 3 for this project

# 0. Read in the data ---- #####################################################
install.packages("triplesmatch")
library(triplesmatch)
load("football_data_prepped_no_outcomes.RData")


# 1. Create 8 strata ---- ######################################################
# Based on propensity score quartile x participation in school government
ps_quart <-  cut_quant(v = logit_ps,  q = c((1:3) / 4),  int = TRUE)
st_factor <-  factor(ps_quart):factor(covs_imputed$schgovt)
st <-  as.integer(st_factor)
names(st) <- names(z)


# 2. Make distance matrix  ---- ################################################

# Discretize Duncan socioeconomic index of father's job into 3 categories
ses3 <-  cut_quant(covs_imputed$bmfoc2u, q = (1:2) / 3)
# Discretize propensity score into 8 categories
logit_ps8 <-  cut_quant(logit_ps, (1:7) / 8)
# Look at covariates of HS rank, Father's education, Mother's education, IQ,
#     propensity score, 3 categories of Duncan SEI of father's job,
#     and 8 categories of propensity score
lcovs <-  cbind(covs_imputed[, c("hsrankq", "bmfaedu", "bmmaedu", "gwiiq_bm")],
                logit_ps, ses3, logit_ps8)
rm(ses3, logit_ps8)
row.names(lcovs) <- row.names(covs_imputed)
# Use those covs for the distance matrices
cost <-  dist_mahal(X = lcovs, st = st, z = z, rank_based = FALSE)


# 3. Run the triples algorithm  ---- ###########################################

m <-  triples(cost, z, st)

# What is the overall objective achieved? 1495.00
m$obj

# What is the lower bound using the full match? 1341.57
m$bound

# How many individuals were matched?
used <- unlist(m$m[, c("itreated", "jcontrol", "kthird")])
sum(duplicated(used)) # Double check that all matched individuals are unique
n <- length(used)
n
# 978

# 4. Make figures and tables of match ---- #####################################

#### Figure 1: Series of three plots for propensity scores ----

setEPS()
postscript("fig1_v4.eps", width = 6, height = 4)

boxplot_matches(m = m$m, y = logit_ps, z = z, yname = "Propensity score")

dev.off()

#### Figure 2: Series of six plots for distributions of propensity score, ----
#    HS rank, Father's education, Mother's education, IQ, and SEI of Father's occupation
setEPS()
postscript("fig2_v3.eps", width = 6.5, height = 5)

par(mfrow=c(2, 3))
par(mar=c(3, 3, 2, 2.5))
boxplot_diffs(m = m$m, y = logit_ps, z = z, yname = "Propensity score")
ycovs <- c("hsrankq", "bmfaedu", "bmmaedu", "gwiiq_bm", "bmfoc2u")
for (cov in ycovs) {
  ycov <- covs_raw[, cov]
  names(ycov) <- row.names(covs_raw)
  boxplot_diffs(m = m$m, y = ycov, z = z, yname = cov_names[cov])
}
par(mar=c(5,  4,  4,  2)  +  0.1)

dev.off()

#### Table 1: Covariate balance  ----

library(xtable)
X<-cbind(logit_ps,ps_quart,covs_raw)
bal <-  make_bal_tab(X = X, z = z, m = m$m,
                   cov_names = c("Propensity Score", "PS-Strata", cov_names))
bal <-  bal[order( - abs(bal[, 4])), ]
rownames(bal)[c(1, 2, 3, 6, 15, 24, 28, 38)] <- paste0("\\textbf{", rownames(bal)[c(1, 2, 3, 6, 15, 24, 28, 38)], "}")
print(xtable(bal, digits=c(0, 2, 2, 2, 2, 2),
             caption = "Covariate balance before/after matching.  
             Means are T = treated, Cb = controls before, 
             Cm = controls after matching and weighting.  
             $\\mathrm{Db = (T-Cb)/s}$ and $\\mathrm{Dm = (T-Cm)/s}$ 
             where s is the pooled st. dev. before matching. 
             \\textbf{Bold} is a covariate used in matching or $|\\mathrm{D}| \\geq 0.2$."),
      caption.placement = "top",
      add.to.row = list(pos = list(0), command = "& $T$ & $Cb$ & $Cmw$ & $Db$ & $Dmw$ \\\\"),
      include.colnames = FALSE,
      sanitize.rownames.function = identity)

#### Tables 2 and 3: Various counts ----

# Make top section of Table 2 with stratum counts
st_factor_collapsed <- factor(st_factor, levels = levels(st_factor), 
                              labels = c("1", "1", "2", "2", levels(st_factor)[5:8]))
tab <- table(z, st_factor_collapsed)[2:1, ]
tab2 <- table(ps_quart)
temp <-  addmargins(tab)
rownames(temp) <-  c("Football", "Control", "Total")
colnames(temp) <-  c("Q1", "Q2",
                     "Q3 - G", "Q3 + G", "Q4 - G", "Q4 + G", "Total")
# Top of table 3 and bottom of table 2
triple_counts <- table(m$m$nOfTreated, m$m$st)
triple_counts[, 2] <- triple_counts[, 2] + triple_counts[, 3]
triple_counts <- triple_counts[, -3]
matched_football <- sapply(list(1:2, 3:4, 5, 6, 7, 8), 
                           function(x) sum(m$m$nOfTreated[m$m$st %in% x]))
matched_control <- colSums(triple_counts) * 3 - matched_football
matched_total <- matched_football + matched_control
temp2 <- rbind(triple_counts, matched_football, matched_control, matched_total)
temp2 <- cbind(temp2, rowSums(temp2))
rownames(temp2) <-  c("1 - 2 triples", "2 - 1 triples",
                      "Matched Football", "Matched Control", "Matched Total")

# Bottom section of Table 3 with flow info
mtk <- apply(tab, 2, function(x) (min(x[1], 2*x[2])))
mck <- tab[2, ]

tab_full <- table(z, st_factor)[2:1, ]
bk <- pmax(0, ceiling((2 * apply(tab_full, 2, function(x) (min(x[1], 2*x[2]))) - 
                         tab_full[2, ]) / 3))
bk[1] <- sum(bk[1:2])
bk[3] <- sum(bk[3:4])
bk <- bk[c(1,3,5:8)]
supply <- mtk - bk
temp3 <- rbind(mtk, bk, supply)
rownames(temp3) <- c("Treated to use $= m^*_{tk}$",
                     "2 - 1 triples = $b_k$", "Supply = $m^*_{tk} - b_k$")
temp3 <- cbind(temp3, rowSums(temp3))
temp3 <- rbind(c("1", "2, 3", "4", "5", "6", "7", ""), temp3)
rownames(temp3)[1] <- "Stratum in the Algorithm"
temp3[2, 7] <- paste0("$T^* =$ ", temp3[2, 7] )
temp3[1, 7] <- NA

# Put the components together for table 2
print(xtable(rbind(temp, temp2[3:5, ]), digits = 0,
             align = "r|r|r|rr|rr|r",
             caption = "Counts for each stratum, which is defined by the propensity 
             score quartile (Q1-Q4) and participation in school government (--G and +G). 
             For the sole purpose of a display that complies with WLS rules on confidentiality, 
             the strata for the first two quartiles have been collapsed across school government 
             to avoid reporting small cell counts.",
             label = "tabPropensity"),
      caption.placement = "top",
      include.colnames = FALSE,
      sanitize.text.function = identity,
      hline.after = c(-1, 0, 2, 3, 5, 6),
      add.to.row = list(pos = list(0, 3),
                        command = c("Stratum & Q1 & Q2 & Q3 $-$ G & Q3 $+$ G & Q4 $-$ G & Q4 $+$ G & Total \\\\ \n",
                                    paste0("Quartile Total & \\multicolumn{1}{|c|}{", tab2[1],
                                           "} & \\multicolumn{1}{c|}{", tab2[2],
                                           "} & \\multicolumn{2}{c|}{", tab2[3],
                                           "}  & \\multicolumn{2}{c|}{", tab2[4],
                                           "}  & ", sum(tab2), " \\\\ \n"))))

# Put the components together for table 3
print(xtable(rbind(temp2[1:2, ], temp3), digits = 0,
             align = "r|r|r|rr|rr|r",
             caption = "Counts of 1-to-2 and 2-to-1 matched sets in the triples design.",
             label = "tabPropensity2"),
      caption.placement = "top",
      include.colnames = FALSE,
      sanitize.text.function = identity,
      hline.after = c(-1, 0, 2, 3, 6),
      add.to.row = list(pos = list(0),
                        command = c("Stratum & Q1 & Q2 & Q3 $-$ G & Q3 $+$ G & Q4 $-$ G & Q4 $+$ G & Total \\\\ \n")))





rm(bal, temp, temp2, temp3, triple_counts, bk, cov,
   matched_control, matched_football, matched_total,
   mck, mtk, supply, tab, tab2, ycov, ycovs)


# 5. Look at marginal outcomes ---- ###############################################

# Before doing the analysis, we can look at marginal outcomes
# Load in outcome data
load("football_data_prepped_w_outcomes.RData")

# Remove treated vector since we can't use it
rm(z)

# Reduce down to only matched individuals
outcomes <- outcomes[used, ]


## Create ranks for matched individuals ----

### Create aberrant ranks based on adjusted TICSm scores alone ----

# How many matched individuals had TICSm below 29?
sum(outcomes$ticsm_adj < 29)
sum(outcomes$ticsm_adj < 29) / n * 100
# 315 aberrant individuals, or 32.2%

# Anyone above 29 receives a rank of 0
# The lower the TICSm score, the higher the aberrant rank
outcomes$ab <- aberrantscoreslong(y = outcomes$ticsm_adj, cutoff = 29, cutoff_dir = "less")
table(outcomes$ab)
# Smallest positive aberrant ranks go to individuals with TICSm of 28
table(outcomes$ticsm_adj[outcomes$ab == min(outcomes$ab[outcomes$ab > 0])])
# Highest aberrant rank of 314.5 is shared by two individuals with TICSm of 7
max(outcomes$ab)
table(outcomes$ticsm_adj[outcomes$ab == max(outcomes$ab)])


### Create coherent ranks based on TICSm and diagnoses ----

outcomes$diagnosis_num <- as.numeric(factor(outcomes$diagnosis,
                                            levels = c("Dementia", "MCI", "Normal"),
                                            labels = 0:2, exclude = NA))
aberrant <- outcomes$ab > 0

# What percent of aberrant diagnoses are available?
sum(!is.na(outcomes$diagnosis_num[aberrant]))
sum(!is.na(outcomes$diagnosis_num[aberrant])) / sum(aberrant) * 100

# Compare TICSm scores of each pair of individuals
# 1 if worse, -1 if better
outcomes$coh <- NA
ticsm_adj_mat_g <- outer(outcomes$ticsm_adj[aberrant], outcomes$ticsm_adj[aberrant], ">")
ticsm_adj_mat_l <- outer(outcomes$ticsm_adj[aberrant], outcomes$ticsm_adj[aberrant], "<")
ticsm_adj_mat <- 1 * ticsm_adj_mat_l - 1 * ticsm_adj_mat_g
sum(ticsm_adj_mat)

# Compare diagnoses of each pair of individuals
# 1 is worse, -1 if better
# 0 if one is missing diagnosis
diagnosis_mat_g <- outer(outcomes$diagnosis_num[aberrant], outcomes$diagnosis_num[aberrant], ">")
diagnosis_mat_l <- outer(outcomes$diagnosis_num[aberrant], outcomes$diagnosis_num[aberrant], "<")
diagnosis_mat <- 1 * diagnosis_mat_l - 1 * diagnosis_mat_g
sum(diagnosis_mat, na.rm = TRUE)
diagnosis_mat_nona <- diagnosis_mat
diagnosis_mat_nona[is.na(diagnosis_mat_nona)] <- 0

# Combine comparisons of ticsm and diagnoses
# Only need strict inequality for one of the two comparisons,
#      can be equal or nonexistant for the other
coherent <- 1 * ((ticsm_adj_mat + diagnosis_mat_nona) >= 1) -
  1 * ((ticsm_adj_mat + diagnosis_mat_nona) <= -1)

# Create scores for each individual by summing all their comparisons
outcomes$coh[aberrant] <- rowSums(coherent)

# Turn these into ranks by shifting such that smallest rank is 1
outcomes$coh[aberrant] <- outcomes$coh[aberrant] -
  min(outcomes$coh[aberrant]) + 1
# Nonaberrant individuals receive coherent ranks of 0
outcomes$coh[!aberrant] <- 0

rm(diagnosis_mat, diagnosis_mat_g, diagnosis_mat_l, diagnosis_mat_nona,
   ticsm_adj_mat, ticsm_adj_mat_g, ticsm_adj_mat_l, coherent)

# How do the coherent ranks compare to aberrant ranks of TICSm, TICSm, and diagnoses?
cor(outcomes$coh[aberrant], outcomes$ab[aberrant])
cor(outcomes$coh[aberrant], outcomes$ticsm_adj[aberrant])
cor(outcomes$coh[aberrant], outcomes$diagnosis_num[aberrant], use = "complete.obs")


## Plots for protocol ----

library(ggplot2)
library(dplyr)
coh_df <- outcomes[aberrant, ] %>% group_by(coh, ticsm_adj, diagnosis)
coh_df <- coh_df %>% summarize(count = n(), .groups = "keep")

## Scatterplot: coherent scores vs TICSm and diagnoses among matched aberrant ----

setEPS()
postscript("coherent_v2.eps", width = 8.0, height = 6.0)

ggplot(coh_df, mapping = aes(y = coh, x = ticsm_adj,
                             col = diagnosis, size = count)) +
  geom_point() + theme_bw() +
  xlab("Adjusted TICSm score") +
  ylab("Coherent score") +
  ggtitle("Coherent score vs components among aberrant for matched individuals")

dev.off()

rm(coh_df)


## Boxplot: coherent scores vs diagnoses among matched aberrant ----

setEPS()
postscript("coherent_boxplot_v2.eps", width = 8.0, height = 6.0)

boxplot(outcomes$coh[aberrant] ~ outcomes$diagnosis[aberrant],
        ylab = "Coherent score", xlab = "Diagnosis",
        main = "Coherent score vs diagnosis among aberrant for matched individuals")

dev.off()



save.image("football_match_081323.RData")
