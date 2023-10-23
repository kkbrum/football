# This script is used for the final analysis (and to create relevant figures and tables)
# This is script 3 of 3 for this project


library(informedSen)
library(sensitivityfull)
library(xtable)

# 0. Load data ---- ############################################################

load("football_data_prepped_w_outcomes.RData")
load("football_match_081323.RData")
install.packages("triplesmatch")
library(triplesmatch)

z <- z[used]

# 1. Plot outcomes vs treatment ---- ###########################################

m12 <- unlist(m$m[m$m$nOfTreated == 1, c("itreated", "jcontrol", "kthird")])
mtype <- as.numeric(rownames(outcomes) %in% m12)
mtype <- factor(mtype, levels = c(1, 0), labels = c("1-2", "2-1"))
zfactor <- factor(z, levels = c(1, 0), labels = c("T", "C"))

setEPS()
postscript("fig4.eps", width = 10.0, height = 6.0)

par(mfrow = c(1,4))
boxplot(outcomes$ticsm_adj ~ zfactor + mtype, 
        ylab="TICSm", xlab="Group", main="(i) TICSm",
        las=1, names = paste(rep(levels(zfactor), 2), 
                             rep(levels(mtype), each = 2), sep = " "))
boxplot(outcomes$ab ~ zfactor + mtype, 
        ylab="Aberrant Ranks", xlab="Group",
        main="(ii) Aberrant TICSm",las=1, 
        names = paste(rep(levels(zfactor), 2), 
                      rep(levels(mtype), each = 2), sep = " "))
boxplot(outcomes$coh ~ zfactor + mtype, 
        ylab="Coherent Ranks", xlab="Group",
        main="(iii) Coherent",las=1,
        names = paste(rep(levels(zfactor), 2), 
                      rep(levels(mtype), each = 2), sep = " "))
boxplot(c(outcomes$coh[z == 1 & mtype == "1-2"], 
          rep(outcomes$coh[z == 1 & mtype == "2-1"], 2)),
        c(outcomes$coh[z == 0 & mtype == "1-2"], 
          rep(outcomes$coh[z == 0 & mtype == "2-1"], 2)), 
        main="(iv) Weighted", 
        ylab="Coherent Ranks",
        names=c("T","C"), xlab="Group"
)

dev.off()

# 2. Tests for effects ---- ####################################################

## TICSm: no effect vs harmful (lower TICSm) ----

y_ticsm <- outcomes$ticsm_adj
names(y_ticsm) <- row.names(outcomes)
data_sen_ticsm <- formattrip(m = m$m, y = y_ticsm)
data_sen_ticsm$zmat <- cbind(data_sen_ticsm$treated1, 
                                !data_sen_ticsm$treated1,
                                !data_sen_ticsm$treated1)

# Do treated people have lower TICSm scores?
sentrip(scores = data_sen_ticsm$ymat, treated1 = data_sen_ticsm$treated1,
        gamma = 1, alternative = "less")
# P-value 0.763

## TICSm: decline of 2 points vs less harmful (higher TICSm) ----
sentrip(scores = data_sen_ticsm$ymat + 2 * data_sen_ticsm$zmat, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1, alternative = "greater")
# P-value 5.63e-12

# With bias of Gamma = 2
sentrip(scores = data_sen_ticsm$ymat + 2 * data_sen_ticsm$zmat, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 2, alternative = "greater")
# P-value 0.0199

# With bias of Gamma = 2.13
sentrip(scores = data_sen_ticsm$ymat + 2 * data_sen_ticsm$zmat, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 2.13, alternative = "greater")
# P-value 0.0515

## Aberrant ranks: no effect vs harm (higher rank) ----

y_ab <- outcomes$ab
names(y_ab) <- row.names(outcomes)
data_sen_ab <- formattrip(m = m$m, y = y_ab)

sentrip(scores = data_sen_ab$ymat, 
        treated1 = data_sen_ab$treated1,
        gamma = 1, alternative = "greater")
# P-value 0.77

## Aberrant ranks: decline of 2 on TICSm scale vs less harmful (lower rank) ----

# Recompute aberrant ranks after accounting for the null effect of -2 for treated
y_ab2 <- aberrantscores(data_sen_ticsm$ymat, cutoff = 29, cutoff_dir = "less",
                        tau = -2, treated1 = data_sen_ticsm$treated1)

sentrip(scores = y_ab2, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1, alternative = "less")
# P-value 7.97e-06

# With bias Gamma = 1.6
sentrip(scores = y_ab2, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1.6, alternative = "less")
# P-value 0.038

# With bias Gamma = 1.64
sentrip(scores = y_ab2, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1.64, alternative = "less")
# P-value 0.0503

## Coherent ranks: no effect vs harm (greater rank) ----

y_coh <- outcomes$coh
names(y_coh) <- row.names(outcomes)
data_sen_coh <- formattrip(m = m$m, y = y_coh)

# Do treated people have lower TICSm scores?
sentrip(scores = data_sen_coh$ymat, treated1 = data_sen_coh$treated1,
        gamma = 1, alternative = "greater")
# P-value 0.713


# 3. Equivalence interval ---- #################################################

# What is the greatest decline that is not rejected when testing the TICSm score?
# Test versus less harm, which is greater TICSm scores
sentrip(scores = data_sen_ticsm$ymat + 0.29 * data_sen_ticsm$zmat, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1, alternative = "greater")
# P-value 0.0514

sentrip(scores = data_sen_ticsm$ymat + 0.3 * data_sen_ticsm$zmat, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1, alternative = "greater")
# P-value 0.048

# We did not reject 0 as too small a harm earlier, 
#   so the one-sided interval is [-0.29, 0]

## At Gamma = 1.25:

sentrip(scores = data_sen_ticsm$ymat, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1.25, alternative = "less")
# P-value 0.99 for testing 0 vs harm

sentrip(scores = data_sen_ticsm$ymat + 0.78 * data_sen_ticsm$zmat, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1.25, alternative = "greater")
# P-value 0.053

sentrip(scores = data_sen_ticsm$ymat + 0.79 * data_sen_ticsm$zmat, 
        treated1 = data_sen_ticsm$treated1,
        gamma = 1.25, alternative = "greater")
# P-value 0.049

# One-sided interval at Gamma = 1.25 is [-0.78, 0]


# 4. Check for bias with unaffected outcomes ---- ##############################

y_periodont <- outcomes$periodont
names(y_periodont) <- row.names(outcomes)
data_mh_periodont <- formattrip(m = m, y = y_periodont, type = "long")
mantelhaen.test(x = table(data_mh_periodont$z, data_mh_periodont$y, data_mh_periodont$mset),
                alternative = "two.sided")
# P-value of 0.26

y_cancer <- outcomes$cancer
names(y_cancer) <- row.names(outcomes)
data_mh_cancer <- formattrip(m = m, y = y_cancer, type = "long")
mantelhaen.test(x = table(data_mh_cancer$z, data_mh_cancer$y, data_mh_cancer$mset),
                alternative = "two.sided")
# P-value of 0.16

# We don't find any bias

