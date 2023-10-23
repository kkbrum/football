# This script is used prepare the data that will be used for matching and for analyses
# This is script 1 of 3 for this project

# 0. Read in the data ---- #####################################################

# Public versions of the WLS data can be downloaded at https://researchers.wls.wisc.edu/data/survey-data/
# This analysis requires the private data, which can be requested by following the instructions at the link above

# This first data set is WLS version 14.01 in wide format, including the private id variable (idpriv)
# The private data is needed in order to access the yearbook activity data, including football playing
# This was obtained 01/30/2023 by Dylan Small from Carol Roan
data_football <- haven::read_dta("Small_14_01_Football_Complete.dta")
data_football <- haven::zap_label(data_football)
data_football <- haven::zap_formats(data_football)
data_football <- haven::zap_labels(data_football)

# This second data set is WLS version 14.02 in wide format, including the private id variable (idpriv)
# This data is needed because it adds the results from the long interviews
# The private id is needed in order to merge this to the first data set
# This was obtained 08/24/2023 by Dylan Small from Carol Roan
data_1402 <- haven::read_dta("Small_14_02_idpriv.dta")
data_1402 <- haven::zap_label(data_1402)
data_1402 <- haven::zap_formats(data_1402)
data_1402 <- haven::zap_labels(data_1402)

## Combine the first and second data sets
# Some variables were updated a bit between the two data sets, keep the newer version of everything available
data_football <- data_football[, (!names(data_football) %in% names(data_1402)) | (names(data_football) == "idpriv")]
data <- merge(data_1402, data_football, by = "idpriv")
rm(data_football, data_1402)


# 1. Eligibility and exclusion criteria ---- ###################################


## 10317 people total ----
nrow(data)

## Restrict to 4991 men (sexrsp == 1) ----
table(data$sexrsp)
data <- data[data$sexrsp == 1, ]

## Yearbook information  ----
# https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=ancillary&module=actv 

### Drop 621 students with no yearbook containing activity info ----

# 0 = no yearbook, 
# 1 = no activity info in yearbook, 
# 2 = no activity info for that student in their yearbook
# 3 = complex yearbook with no activity info for that student
# 4 = student not in yearbook
table(data$activ1[data$activ1 %in% 0:4])
sum(data$activ1 %in% c(0, 1, 4))
data <- data[!data$activ1 %in% c(0, 1, 4), ]

### Drop 397 students with complex yearbook ----
sum(data$cmplxschl == 1)
data <- data[data$cmplxschl != 1, ]

### Left with 3973 men with adequate yearbook info ----

## Sports participation ----

data$baseball <- 0
data$basketball <- 0
data$crosscountry <- 0
data$curling <- 0
data$football <- 0
data$lacrosse <- 0 
data$soccer <- 0
data$track <- 0
data$volleyball <- 0
data$wrestling <- 0
data$swimming <- 0
data$hockey <- 0
data$gymnastics <- 0
data$tennis <- 0
data$othersport <- 0

sports <- c(baseball = 1001, basketball = 1002, crosscountry = 1003,
            curling = 1004, football = 1005, lacrosse = 1006, 
            soccer = 1007, track = 1008, volleyball = 1009, wrestling = 1010,
            swimming = 1011, hockey = 1012, gymnastics = 1013, tennis = 1014,
            othersport = 1015)

# Look across the 27 activity variables
for (i in 1:27) {
  # Look across all possible sports
  for (j in 1:length(sports)) {
    # Mark down participation
    idx <- which(data[, paste0("activ", i)] == sports[j])
    data[idx, names(sports[j])] <- 1
  }
}

# What sports did the males with proper yearbook info play?
colSums(data[, names(sports)])

### 2627 people did not play football ----
sum(data$football == 0)

### Exclude 65 people who play collision sports and not football ----
data$collision <- 0
data$collision[rowSums(data[, c("wrestling", "hockey", "lacrosse", "soccer")]) >= 1] <- 1
sum(data$collision[data$football == 0])
data <- data[data$collision == 0 | data$football == 1, ]

### Now have 3908 eligible subjects: 1346 football and 2562 non-football  ----
table(data$football)

# What other sports did these groups play?
colSums(data[data$football == 1, names(sports)])
colSums(data[data$football == 0, names(sports)])

rm(i, idx, j, sports)



# 2. Outcomes ---- #############################################################

## Check availability of TICS-M measured in the 2020 wave ----
# https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=iliadgradw1&module=q1tics

# Age at this wave 
table(data$q1a003re) # age at short interview
range(data$q1a003re + data$q1a941re/12, na.rm = TRUE) # approximate age at long interview

# q1i915re is the raw TICSm score
data$ticsm_raw <- data$q1i915re

# Completion of short interview (which yields TICSme)
# 1: complete, 2: refused, 3: not found, 4: deceased, 5: respondent away/unavailable, 6: unable
table(data$stat20short)

### Left with 1156 men with primary outcome data; 729 control and 427 treated ----
include <- !is.na(data$ticsm_raw)
sum(include)
table(data$football[include])


## Check availability of negative outcomes ----


### Periodontitis at 81 wave: 1 refused, 5 unknown ----

# https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=iliadgradw1&module=q1health
# ``Has a doctor ever told Participant they have periodontitis or gum disease?"
# -3 is refused, -1 is don't know, 1 is yes, 2 is no.
data$periodont <- data$q1x3991re

# Missingness is the same as for TICS-M since it is in the same wave
sum(!is.na(data$periodont)) 
table(data$periodont)

# Count unknown as "no" due to wording of the question
data$periodont[data$periodont == -1] <- 2

#### Exclude the person who refused neg outcome ----
include <- include & data$periodont != -3

# Recode to 0/1
data$periodont[data$periodont == -3] <- NA
data$periodont[data$periodont == 2] <- 0




### Cancer at 81 wave: 1 refused ----
# https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=iliadgradw1&module=q1health
# ``Has a doctor ever told Participant they have cancer or a malignant tumor, not including minor skin cancers?"
# -3 is refused, 1 is yes, 2 is no.
data$cancer <- data$q1x348re

# Missingness is the same
sum(!is.na(data$cancer)) 
table(data$cancer)

# The person that refused periodontitis also refused cancer, so we already drop them
table(data$cancer[include])

## Left with 1155 men; 426 football and 729 not
sum(include)
table(data$football[include])

# Recode to 0/1
data$cancer[data$cancer == -3] <- NA
data$cancer[data$cancer == 2] <- 0



## Long interview data ----
# https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=iliadgradw1&module=q1beg

### Create adjusted TICSm score ----
data$ticsm_adj <- data$ticsm_raw - 2 * data$q1i914re - 1 * data$q1i913re
table(data$ticsm_adj)
sum(data$ticsm_adj[include] < 29)

# Long interview participation
# 1 is completed, 2 is not-eligible, 3 is refused/noncontact,
# 4 is died before completing, 5 is physically/mentally unable
table(data$stat20long[include & data$ticsm_adj < 29])

# Dementia questionnaire participation
# 1 is completed, 2 is completed for dead people, and 3 is not completed
table(data$stat20dq[include & data$ticsm_adj < 29], 
      data$stat20long[include & data$ticsm_adj < 29])

### Create diagnosis ----

# Use the positive predictive value diagnosis
# 1 is > TICSm cutoff
# 2 and 3 are normal cognition via consensus vs proxy
# 4 and 5 are MCI with and without Alzheimer's (AD)
# 6, 7, 8, and 9 are dementia via consensus with AD, 
#           consensus without AD, proxy with AD, proxy without AD
# 11 is died before or refused long interview
table(data$q1a956re[include & data$ticsm_adj < 29])

# 2 people were physically/mentally unable to complete long interview and did not get a proxy
sum(is.na(data$q1a956re[include & data$ticsm_adj < 29]))
data[is.na(data$q1a956re) & include & data$ticsm_adj < 29, 
     c("ticsm_adj", "stat20long", "stat20dq")]
# Give these people a score of 12
data$q1a956re[is.na(data$q1a956re) & include & data$ticsm_adj < 29] <- 12

data$q1a956re[data$q1a956re == 1] <- NA
table(data$q1a956re[include & data$ticsm_adj < 29])

data$diagnosis <- factor(data$q1a956re, levels = c(2:9, 11, 12), 
                         labels = c(rep("Normal", 2), rep("MCI", 2),
                                    rep("Dementia", 4), rep("Died/refused/unable", 2)))

table(data$diagnosis[include & data$ticsm_adj < 29])

# Who has a diagnosis but shouldn't (TICSm >= 29)? 
accidental_diagnosis <- include & !is.na(data$diagnosis) & data$ticsm_adj >= 29
data[accidental_diagnosis, c("ticsm_adj", "stat20long", "stat20dq", "diagnosis")]
# stat20long = 1 means the long survey was conducted, 2 means not eligible
# stat20dq = 3 means not completed, 2 means proxy completed because participant died

# Take out these extra diagnoses
data$diagnosis[accidental_diagnosis] <- NA
rm(accidental_diagnosis)


# 3. Define covariates ---- ####################################################


cov_names <- c("spocasp3" = "Duncan socioeconomic index of aspired job", 
               "hsrankq" = "High school rank", 
               "bmpin1" = "Parental income", 
               "bmfaedu" = "Father education level", 
               "bmmaedu" = "Mother education level", 
               "bmfoc2u" = "Duncan socioeconomic index of father's job", 
               "zpedyr" = "Planned future education",
               "gwiiq_bm" = "IQ", 
               "tchevl" = "Outstanding student", 
               "musperf" = "Participated in musical ensemble", 
               "spchperf" = "Participated in debate",
               "schgovt" = "Participated in school govt", 
               "schpubs" = "Participated in school pubs", 
               "rlur57" = "Rural: father was a farmer",
               "milit" = "Planned to serve in military", 
               "cath" = "Catholic high school", 
               "bklvpr" = "Lived with both parents", 
               "wrmo57" = "Mother working", 
               "tchencq" = "Teachers encouraged college", 
               "parencq" = "Parents encouraged college", 
               "zfrplc" = "Friends planned on college", 
               "tchcntq1" = "Discussed plans with teachers: Not at all", 
               "tchcntq2" = "Discussed plans with teachers: Sometimes", 
               "tchcntq3" = "Discussed plans with teachers: Very much", 
               "parcntq1" = "Discussed plans with parents: Not at all", 
               "parcntq2" = "Discussed plans with parents: Sometimes", 
               "parcntq3" = "Discussed plans with parents: Very much", 
               "sesp571" = "Family wealth: Considerably below", 
               "sesp572" = "Family wealth: Somewhat below", 
               "sesp573" = "Family wealth: Average", 
               "sesp574" = "Family wealth: Somewhat above", 
               "sesp575" = "Family wealth: Considerably above", 
               "parsup1" = "Parent support for college: cannot",
               "parsup2" = "Parent support for college: can with sacrifice",
               "parsup3" = "Parent support for college: can easily",
               "ixa04rec" = "Family history of early stroke",
               "ixa09rec" = "Family history of alzheimer's",
               "ixc02rer" = "Childhood asthma",
               "ixc07rer" = "Childhood polio")

# Recode musperf, spchperf, schgovt, schpubs to no/yes
data$musperf <- as.numeric(data$musperf >= 1)
data$spchperf <- as.numeric(data$spchperf >= 1)
data$schgovt <- as.numeric(data$schgovt >= 1)
data$schpubs <- as.numeric(data$schpubs >= 1)

# Recode plns58q to military no/yes
# See https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=wls5764&module=aasp for definitions
data$milit <- as.numeric(data$plns58q %in% c(4,11,17,22,27,28,30))

# Recode hsdm57 to catholic no/yes
data$cath <- as.numeric(data$hsdm57 == 1)

# Recode rural/farmer to no/yes
# See https://www.ssc.wisc.edu/wlsresearch/documentation/waves/?wave=wls5764&module=asum for exact definition
data$rlur57 <- as.numeric(data$rlur57 == 1)

# zpedyr value of -2 value of "inappropriate" means they did not plan on any more schooling
data$zpedyr[data$zpedyr == -2] <- 0

# Recode continuous variables with missingness indicators
covs_cont_w_missing <- c("spocasp3", "hsrankq", "bmpin1", "bmfaedu", "bmmaedu",
                         "bmfoc2u", "zpedyr")
for (cov in covs_cont_w_missing) {
  data[, paste0(cov, "missing")] <- as.numeric(data[, cov] < 0 | is.na(data[, cov]))
  data[data[, paste0(cov, "missing")] == 1, cov] <- NA
}

# Recode some categorical variables to missing/no/yes
covs_cat_w_missing <- c("bklvpr", "wrmo57", "tchencq", "parencq", "zfrplc",
                        "ixa04rec", "ixa09rec", "ixc02rer", "ixc07rer")
for (cov in covs_cat_w_missing) {
  data[, paste0(cov, "missing")] <- as.numeric(data[, cov] < 0 | is.na(data[, cov]))
  data[, cov] <- as.numeric(data[, cov] == 1)
  data[data[, paste0(cov, "missing")] == 1, cov] <- NA
}

# Recode ordinal variables as dummy variables and add missing indicator (although it will be colinear)
covs_ordinal_w_missing <- c("tchcntq", "parcntq", "sesp57", "parsup")
for (cov in covs_ordinal_w_missing) {
  if (cov == "sesp57") { 
    ncat <- 5 
  } else {
    ncat <- 3
  }
  for (i in 1:ncat) {
    data[, paste0(cov, i)] <- as.numeric(data[, cov] == i)
    data[data[, cov] == -3, paste0(cov, i)] <- NA
  }
  data[, paste0(cov, 1, "missing")] <- as.numeric(data[, cov] == -3)
  data[data[, cov] == -3, cov] <- NA
}
rm(cov, i, ncat)



# 4. Create final sample ---- ##################################################

## Separate out deceased/unavail data ----
deceased_data <- data[!is.na(data$stat20short) & data$stat20short == 4, ]
# People who are invited for this wave but are missing primary/neg outcome and are not dead
unavail_data <- data[!is.na(data$stat20short) & !include & data$stat20short != 4, ]

## Sample includes eligible people who have primary and negative outcomes available ----
sample_data <- data[include, ]
rm(data, include)


# 5. Covariate missingness ---- ################################################

## Table of covariate missingness ----

missingness_vars <- c(covs_cont_w_missing, covs_cat_w_missing, 
                      paste0(covs_ordinal_w_missing, 1))
missingness_cov_names <- cov_names[missingness_vars]
missingness_cov_names <- sapply(missingness_cov_names, function(x) {
  strsplit(x, split = ":")[[1]][1]})
missingness_cov_names <- paste0(missingness_cov_names, "- missing")
names(missingness_cov_names) <- paste0(missingness_vars, "missing")

missingness <- colSums(is.na(sample_data[, names(cov_names)]))
names(missingness) <- cov_names
xtable::xtable(matrix(missingness, ncol = 1, dimnames = list(cov_names)), digits = 0)

## Include indicators for >= 35 missing values ----
missing_to_include <- missingness_cov_names[colSums(sample_data[, names(missingness_cov_names)]) >= 35]
missing_to_exclude <- which(names(missing_to_include) %in% c("ixa09recmissing", "parsup1missing"))
missing_to_include <- missing_to_include[-missing_to_exclude]
rm(missing_to_exclude)

## Impute missing covariate values ----
cov_names <- c(cov_names, missing_to_include)
covs_raw <- sample_data[, names(cov_names)]
covs_imputed <- covs_raw
colSums(is.na(covs_imputed))

### For binary and categorical variables, impute 0's for the missing values
for (cov in c(paste0("tchcntq", 1:3), paste0("parcntq", 1:3),
              paste0("sesp57", 1:5), paste0("parsup", 1:3),
              covs_cat_w_missing)) {
  covs_imputed[is.na(covs_imputed[, cov]), cov] <- 0
}

### For continuous variables, impute the mean for the missing values
for (cov in covs_cont_w_missing) {
  covs_imputed[is.na(covs_imputed[, cov]), cov] <- mean(covs_imputed[, cov], na.rm = TRUE)
}

colSums(is.na(covs_imputed))
rm(cov, covs_cat_w_missing, covs_cont_w_missing, covs_ordinal_w_missing,
   missing_to_include, missingness, missingness_cov_names, missingness_vars)


# 6. Propensity scores ---- ####################################################
z <- sample_data$football
names(z) <- row.names(sample_data)

logit_model <- glm(z ~ ., data = covs_imputed, family = binomial())
logit_ps <- predict(logit_model, type = "response")
rm(logit_model)



# 7. Marginal TICSm/ diagnosis plots ---- ######################################

setEPS()
postscript("marginal_outcome_hist.eps", width = 8.0, height = 6.0)

hist(sample_data$ticsm_adj, breaks = 30,
     main = "Marginal distribution of adjusted TICSm overall",
     xlab = "Adjusted TICSm score")
abline(v = 29, lwd = 2)

dev.off()

setEPS()
postscript("marginal_outcome_boxplot.eps", width = 8.0, height = 6.0)

boxplot(sample_data$ticsm_adj ~ sample_data$diagnosis,
        main = "Marginal distribution of TICSm among aberrant",
        ylab = "Adjusted TICSm score",
        xlab = "Diagnosis")

dev.off()


# 8. Clean up and save ---- ####################################################
outcome_names <- c("ticsm_raw", "ticsm_adj", "diagnosis", 
                   "periodont", "cancer")
outcomes <- sample_data[, outcome_names]

save.image("football_data_prepped_w_outcomes.RData")

sample_data <- sample_data[, ! names(sample_data) %in% outcome_names]
rm(outcomes, outcome_names)

save.image("football_data_prepped_no_outcomes.RData")

