# R code for replicating the dataset and the analyses in the paper #
# Cannabidiol for psychotic disorders: a meta-analysis of randomized controlled trials #
# R code by Mattia Marchi (mattiamarchimd@gmail.com) 
# August 2, 2025

###------------------------------------------------------------------------------------------------------
###---------------------Meta Analysis - Cannabidiol (CBD) in psychotic disorders-------------------------
###------------------------------------------------------------------------------------------------------
#Load required packages
library(meta)
library(metafor)
library(tidyverse)
library(metaforest)
library(brms)
library(bayesplot)

#Import data
df <- structure(list(ID = structure(c(1L, 2L, 4L, 5L, 6L, 7L, 8L, 3L), levels = c("Boggs et al 2018", "Chesney et al 2025", "Jazz Pharmaceuticals, 2023", "Köck et al 2021", "Leweke et al 2012", "McGuire et al 2018", "O’Neill et al 2021", "Ranganathan 2022"), class = "factor"),
                     Trial.ID.es.NCT54268 = structure(c(1L, 6L, 7L, 2L, 3L, NA, 4L, 5L), levels = c("NCT00588731", "NCT00628290", "NCT02006628", "NCT02504151", "NCT04421456", "NCT04605393", "NCT04700930"), class = "factor"),
                     Status = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L), levels = c("Published", "Unpublished"), class = c("ordered", "factor")), Year = c(2018L, 2025L, 2021L, 2012L, 2018L, 2021L, 2022L, 2023L),
                     Country = structure(c(5L, 3L, 2L, 1L, 4L, 3L, 5L, 6L), levels = c("Germany", "Switzerland", "UK", "UK, Romania, Poland", "USA", "USA, Poland, Serbia, Spain"), class = "factor"), 
                     Study.design = structure(c(2L, 1L, 2L, 2L, 2L, 1L, 1L, 2L), levels = c("Crossover RCT", "Parallel RCT"), class = "factor"), Setting = structure(c(3L, 3L, 1L, 1L, 2L, 3L, 2L, 2L), levels = c("Inpatient", "NR", "Outpatient"), class = "factor"),
                     Response = structure(c(2L, 1L, 2L, 2L, 2L, 1L, 2L, 3L), levels = c("Acute", "Early", "Long"), class = c("ordered", "factor")), Followup = structure(c(6L, 2L, 4L, 4L, 6L, 3L, 5L, 1L), levels = c("12 weeks", "20 minutes", "270 minutes", "28 days", "4 weeks", "6 weeks"), class = "factor"), 
                     Diagnosis = structure(c(4L, 4L, 2L, 3L, 2L, 1L, 1L, 3L), levels = c("Psychosis early stage", "Psychotic disorders", "SZ", "SZ/SCHAFF"), class = "factor"), Psychosis.stage = structure(c(4L, 3L, 2L, 2L, 3L, 1L, 1L, 3L), levels = c("Early stage", "Acute", "Stable", "Chronic"), class = c("ordered", "factor")),
                     T_dose.mg. = c(600L, 1000L, 200L, 800L, 1000L, 600L, 800L, 1000L), T = structure(c(1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L), levels = c("Cannabidiol", "Cannabidiol cigarettes"), class = "factor"), C = structure(c(2L, 2L, 3L, 1L, 2L,2L, 2L, 2L), levels = c("Amisupride (up to 800 mg)", "PBO", "Tobacco cigarettes"), class = "factor"),
                     Current.treatment = structure(c(4L, 5L, 1L, 2L, 7L, 3L, 3L, 6L), levels = c("AP", "AP suspended 3 days before", "NR", "stable dose AP (1 month)", "stable dose AP (1 month) + THC exposure before treatment", "stable dose AP (8 weeks)", "stable dose AP (max 1 AP, stable for 1 month)"), class = "factor"),
                     Cannabis.Users = structure(c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_, NA_integer_), levels = c("Excluded", "Included"), class = c("ordered", "factor")), Mean.age = c(47.4, 39.7, 35.1, 30.1, 40.8, 27.7, 26.1, 38.4),
                     X.Female = c(30.6, 7, 29, 17.9, 42, 33.3, 20, 26), Outcome1_psychosis.severity = structure(c(NA, NA, 1L, 1L, 1L, 1L, 1L, 1L), levels = "PANSS", class = "factor"),
                     NT1 = c(NA, NA, 16L, 21L, 42L, 15L, 10L, 19L), Mean1_T = c(NA, NA, 65.73, 30.5, 68.1, 41.53, 66.9, -10.69), SD1_T = c(NA, NA, 15.85, 16.4, 14.79, 11, 7.8, 7.63), NC1 = c(NA, NA, 15L, 21L, 44L, 15L, 10L, 19L), Mean1_C = c(NA, NA, 69.17, 30.1, 71.9, 44.6, 69, -8.74), SD1_C = c(NA, NA, 16.96, 24.7, 15.49, 18.07, 6.32, 7.76),
                     Outcome2_cognitività = structure(c(3L, 2L, NA, NA, 1L, NA, NA, NA), levels = c("BACS", "HVLT-R delayed recall", "MATRICS"), class = "factor"), NT2 = c(18L, 30L, NA, NA, 42L, NA, NA, NA), Mean2_T = c(26.4, 3.5, NA, NA, 32.21, NA, NA, NA), SD2_T = c(12.2, 2.79, NA, NA, 6.042, NA, NA, NA),
                     NC2 = c(18L, 30L, NA, NA, 44L, NA, NA, NA), Mean2_C = c(25.4, 4.8, NA, NA, 32.91, NA, NA, NA), SD2_C = c(12.5, 2.93, NA, NA, 7.158, NA, NA, NA), Outcome3_negative.symptoms = structure(c(NA, 1L, NA, NA, 2L, NA, NA, 1L), levels = c("PANSS-N", "SANS"), class = "factor"), NT3 = c(NA, 30L, NA, NA, 42L, NA, NA, 19L),
                     Mean3_T = c(NA, 17.4, NA, NA, 43.6, NA, NA, -1.57), SD3_T = c(NA, 6.43, NA, NA, 16.54, NA, NA, 2.44), NC3 = c(NA, 30L, NA, NA, 44L, NA, NA, 19L), Mean3_C = c(NA, 16.6, NA, NA, 48.4, NA, NA, -2.33), SD3_C = c(NA, 6.43, NA, NA, 15.75, NA, NA, 2.53), N_death_T = c(0L, NA, 1L, NA, 0L, NA, 0L, 0L), N_death_C = c(0L, NA, 0L, NA, 0L, NA, 0L, 0L),
                     N_dropT = c(2L, NA, 5L, 3L, 1L, NA, 1L, 2L), N_dropC = c(2L, NA, 5L, 1L, 1L, NA, 0L, 1L), NT4 = c(18L, NA, 16L, 21L, 42L, NA, 10L, 24L), NC4 = c(18L, NA, 15L, 21L, 44L, NA, 10L, 26L), N_AE_T = c(4L, NA, NA, NA, 3L, NA, NA, 3L), N_AE_C = c(1L, NA, NA, NA, 7L, NA, NA, 3L), NT5 = c(18L, NA, NA, NA, 42L, NA, NA, 24L), NC5 = c(18L, NA, NA, NA, 44L, NA, NA, 26L)),
                row.names = c(NA, -8L), class = "data.frame")

#Overall description of the database
summary(df)

###----------------------------------------------------------------------------------------
###--------------------------1. Psychosis severity-----------------------------------------
###----------------------------------------------------------------------------------------
ps <- metacont(n.e = NT1, mean.e = Mean1_T, sd.e = SD1_T,
               n.c = NC1, mean.c = Mean1_C, sd.c = SD1_C,
               data = df, studlab = ID, sm = "SMD")
ps
forest(ps, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "CBD", label.c = "Controls",
       label.left = "Favours CBD", label.right = "Favours Controls", allstudies = F)
#Subgroup by Status - after replace it with also "Cannabis.Users", "Psychosis.stage", "Response"
ps_sg <- metacont(n.e = NT1, mean.e = Mean1_T, sd.e = SD1_T,
                  n.c = NC1, mean.c = Mean1_C, sd.c = SD1_C,
                  data = df, studlab = ID, sm = "SMD",
                  byvar = Status)
ps_sg
forest(ps_sg, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       overall = T, overall.hetstat = TRUE, print.subgroup.labels = TRUE,
       label.e = "CBD", label.c = "Controls",
       label.left = "Favours CBD", label.right = "Favours Controls", allstudies = F)
#Sequence trial
df <- df[order(df$Year), ]
ps <- metacont(n.e = NT1, mean.e = Mean1_T, sd.e = SD1_T,
               n.c = NC1, mean.c = Mean1_C, sd.c = SD1_C,
               data = df, studlab = ID, sm = "SMD")
cum <- metacum(ps)
forest(cum, sortvar = df$Year, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       overall = T, overall.hetstat = TRUE, print.subgroup.labels = TRUE,
       label.e = "CBD", label.c = "Controls",
       label.left = "Favours CBD", label.right = "Favours Controls", allstudies = F)

###----------------------------------------------------------------------------------------
###--------------------------------2. Cognition--------------------------------------------
###----------------------------------------------------------------------------------------
co <- metacont(n.e = NT2, mean.e = Mean2_T, sd.e = SD2_T,
               n.c = NC2, mean.c = Mean2_C, sd.c = SD2_C,
               data = df, studlab = ID, sm = "SMD")
co
forest(co, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "CBD", label.c = "Controls",
       label.left = "Favours CBD", label.right = "Favours Controls", allstudies = F)

###----------------------------------------------------------------------------------------
###------------------------------3. Negative symptoms--------------------------------------
###----------------------------------------------------------------------------------------
sn <- metacont(n.e = NT3, mean.e = Mean3_T, sd.e = SD3_T,
               n.c = NC3, mean.c = Mean3_C, sd.c = SD3_C,
               data = df, studlab = ID, sm = "SMD")
sn
forest(sn, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "CBD", label.c = "Controls",
       label.left = "Favours CBD", label.right = "Favours Controls", allstudies = F)

###----------------------------------------------------------------------------------------
###---------------------------4. Safety and Tolerability-----------------------------------
###----------------------------------------------------------------------------------------
#-----Dropouts
do <- metabin(event.e = N_dropT, n.e = NT4,
               event.c = N_dropC, n.c = NC4,
               data = df, studlab = ID, sm = "OR")
do
forest(do, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "CBD", label.c = "Controls",
       label.left = "Favours CBD", label.right = "Favours Controls", allstudies = F)
#-----Side effects
se <- metabin(event.e = N_AE_T, n.e = NT5,
              event.c = N_AE_C, n.c = NC5,
              data = df, studlab = ID, sm = "OR")
se
forest(se, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "CBD", label.c = "Controls",
       label.left = "Favours CBD", label.right = "Favours Controls", allstudies = F)