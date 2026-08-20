---
title: "CIBBRiNA - D6.2: Bycatch estimation toolkit"
author: "WP6 - Morgane Pommier, Ailbhe Kavanagh and David Lusseau"
date: "2026-08-20"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_float: true
    toc_collapsed: true
    toc_depth: 3 
    css: style.css
---

This document outlines a typical application of the CIBBRiNA bycatch estimation toolkit, using an example dataset. The associated report "D6.2 Bycatch assessment toolbox: Methods to generate bias-corrected bycatch rates and total bycatch estimates" is available on the CIBBRiNA website at: <https://cibbrina.eu/results/> The toolkit presented here is largely based on BEAM, the Bycatch Estimation and Assessment Matrix developed by ICES working group on bycatch of ETP species (WGBYC, <https://doi.org/10.17895/ices.pub.32922134>).

The original idea of BEAM was developed by Marjorie Lyssikatos and David Lusseau. The code was developed by David Lusseau, André Moan, Paula Gutiérrez Muñoz, Morgane Pommier, Henrik Pärn, Torbjörn Säterberg, Kim Magnus Bærum and Amaia Astarloa Diaz. WGBYC ToR C members contributed significantly to non-coding development and beta testing, particularly Marjorie Lyssikatos, Ailbhe Kavanagh, Caterina Fortuna, Guðjón Már Sigurðsson, Katja Ringdahl, Sara Königson, Ruth Fernandez, Allen Kingston, Lotte Kindt-Larsen and Carlos Pinto. The efforts to get BEAM to ICES Transparent Assessment Framework (TAF) are led by André Moan. The original BEAM code, used for ICES annual assessment by WGBYC and adapted here for CIBBRiNA, is available at: <https://github.com/dlusseau/BEAM>.

Further developments under CIBBRiNA WP6 were led by Morgane Pommier, Ailbhe Kavanagh and David Lusseau.

# Set-up

## R Libraries


``` r
library(data.table) # data wrangling
library(reshape2) # data wrangling
library(dplyr) # data wrangling
library(glmmTMB) # model fitting
library(metafor) # testing heterogeneity 
library(ggeffects) # BPUE estimation
library(emmeans) # BPUE estimation, don't go above emmeans  #1.11.0 !
library(this.path) # for auto-detecting the directory 


library(doSNOW) # parallel computation
library(doParallel) # parallel computation
library(progress) # progress bars
library(stringdist) # data checking and cleaning (string manipulation)
library(stringi) # data checking and cleaning (string manipulation)
library(stringr) # data checking and cleaning (string manipulation)
library(tcltk)
library(ggplot2) # data visualisation
library(scales) # data visualisation
```

## Working environment and parallel computing

The toolkit is provided as a collection of standalone functions stored in the working directory. These functions can be sourced using the commands below, making them available to use within the R environment.

On large datasets, the bycatch estimation process can be computationally intensive and time-consuming. Using parallel computing can substantially reduce running time. The code below includes the necessary steps to configure and enable parallel processing. However, this remains optional, and all functions can also be run in serial mode.


``` r
setwd(this.path::this.dir()) #Working directory to location of this markdown file

#Load the functions:
source("calc_bpue_CIBBRiNA.R") #Estimate bycatch rate
source("zeros_assessment_CIBBRiNA.R") #Assess whether  a 0 bycatch rate can be accepted, given the monitoring coverage
source("calc_total_CIBBRiNA.R") #Estimate total bycatch
source("calc_partial_CIBBRiNA.R") #When total bycatch cannot not be estimated, estimate partial bycatch
source("reliability_estimation_CIBBRiNA.R") #Assess reliability of the estimates

#Set-up parallel computing to speed things up

detectCores()
ncores<-12 #This number can be changed based on the output of detectCores()

cl <- makeCluster(ncores)
registerDoSNOW(cl)
```

# Input data

To implement the BEAM pipeline to estimate bycatch rate and total bycatch, two essential datasets are required:

-   A dataset of monitoring effort, with associated records of bycatch

-   A dataset for total fshing effort, for the fishery in which bycatch needs to be estimated.

Both datasets can be aggregated to some extent using grouping variables (e.g. effort and bycatch records summed by vessel size, subregion, gear type, season...).

## Monitoring & Byatch data

The `bycatch_data` object contains information about monitoring effort and associated observed bycatch events.

In this example, bycatch estimates are of interest for 8 species. We consider two region of interest, north and south, and three types of gear (gillnets, trawlers, longlines). The large regions of interest can be broken down into smaller area, the fishery takes place year-round, but data includes information about season. Finally the fishery involves both large and small vessels.

The key columns are related to the effort, and counts of incidents/individuals bycaught. And then any associated grouping variables to be used in the assessment, either to provide distinct estimates, or investigate and correct for subgroup variability in bycatch patterns.

Within these strata, data has been aggregated to finer level, capturing variable that may influence bycatch rate and cause heterogeneity. In the present example, those variables are area (within a region), season and vessel_size (large vs. small). Each row in the input data therefore represent a unique combination of species x region x gear x season x area x vessel size, with the associated effort and bycatch records.

The data is stored in the `/data` folder within the main working directory, and can be imported using the commands below. The BEAM library relies extensively on the data.table package and its syntax for data manipulation within functions (Barrett et al. 2025). For this reason, it is recommended to read the data directly using `fread()`. If importing data using another method (e.g. `read_csv()`), the data.frame object should be converted to a data.table format using `setDT()` before using any of the toolkit functions. It is also good practice to ensure any grouping variable is converted to a factor, so they are treated accordingly in subsequent analysis.


``` r
bycatch_data <- fread("data/bycatch_simulated_monitoring_records.csv")
bycatch_data[, names(.SD) := lapply(.SD, as.factor), .SDcols = is.character]
```

Looking at the structure of the data, we can see assess the aggregation levels, and identify the groups to be used in the assessment.

Here for instance, the following variables have been used to define groups over which effort (**days**), and number if individuals bycaught (**n_ind**) have been aggregated:

-   **species** (10 levels)

-   **region** (2 levels)

-   **gear** (3 levels)

-   **area** (5 levels, nested in each region)

-   **season** (4 levels)

-   **vessel size** (2 levels)

In this specific example, species, region and gear will be used to define the **unit of assessment** `analysis_resolution` (i.e. one bycatch estimate is required per each species x region x gear combination) while area, season and vessel size will be used as candidate co-variates (e.g. random effects) to explore **variability** of bycatch patterns within the fleet, and correct overall estimates accordingly.


``` r
str(bycatch_data)
```

```
## Classes 'data.table' and 'data.frame':	2400 obs. of  8 variables:
##  $ species    : Factor w/ 10 levels "A","B","C","D",..: 1 2 3 4 5 6 7 8 9 10 ...
##  $ region     : Factor w/ 2 levels "North","South": 1 1 1 1 1 1 1 1 1 1 ...
##  $ gear       : Factor w/ 3 levels "Gillnets","Longlines",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ season     : Factor w/ 4 levels "fall","spring",..: 4 4 4 4 4 4 4 4 4 4 ...
##  $ vessel_size: Factor w/ 2 levels "large","small": 2 2 2 2 2 2 2 2 2 2 ...
##  $ area       : Factor w/ 5 levels "a","b","c","d",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ days       : num  75 89.1 90.8 53.2 91.9 ...
##  $ n_ind      : int  15 2 1 3 7 0 0 1 0 0 ...
##  - attr(*, ".internal.selfref")=<externalptr>
```

``` r
head(bycatch_data)
```

```
##    species region     gear season vessel_size   area     days n_ind
##     <fctr> <fctr>   <fctr> <fctr>      <fctr> <fctr>    <num> <int>
## 1:       A  North Gillnets winter       small      a 74.95034    15
## 2:       B  North Gillnets winter       small      a 89.11956     2
## 3:       C  North Gillnets winter       small      a 90.77794     1
## 4:       D  North Gillnets winter       small      a 53.23052     3
## 5:       E  North Gillnets winter       small      a 91.91001     7
## 6:       F  North Gillnets winter       small      a 80.65606     0
```

## Fishing data

The `fishing_data` object contains information about the global (e.g. fleet level) fishing effort taking place, for which total bycatch is to be estimated. The data must be aggregated at the same level than the monitoring data, for instance here: region x gear x season x area x vessel size, each row associated with a value of corresponding effort (days).


``` r
fishing_data <- fread("data/total_simulated_fishing_effort.csv")
fishing_data[, names(.SD) := lapply(.SD, as.factor), .SDcols = is.character]

str(fishing_data)
```

```
## Classes 'data.table' and 'data.frame':	2400 obs. of  6 variables:
##  $ region     : Factor w/ 2 levels "North","South": 1 1 1 1 1 1 1 1 1 1 ...
##  $ gear       : Factor w/ 3 levels "Gillnets","Longlines",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ season     : Factor w/ 4 levels "fall","spring",..: 4 4 4 4 4 4 4 4 4 4 ...
##  $ vessel_size: Factor w/ 2 levels "large","small": 2 2 2 2 2 2 2 2 2 2 ...
##  $ area       : Factor w/ 5 levels "a","b","c","d",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ days       : num  65661 47066 43110 49160 42585 ...
##  - attr(*, ".internal.selfref")=<externalptr>
```

``` r
head(fishing_data)
```

```
##    region     gear season vessel_size   area     days
##    <fctr>   <fctr> <fctr>      <fctr> <fctr>    <num>
## 1:  North Gillnets winter       small      a 65661.33
## 2:  North Gillnets winter       small      a 47065.62
## 3:  North Gillnets winter       small      a 43110.07
## 4:  North Gillnets winter       small      a 49159.93
## 5:  North Gillnets winter       small      a 42585.04
## 6:  North Gillnets winter       small      a 31650.15
```

# Bycatch Rate

## BPUE estimation

The `calc_bpue()` function estimates bycatch rates from an input monitoring dataset, assesses underlying heterogeneity in rates across subgroups (if any), compares all possible models and identifies the best one. It also returns the population-level value of BPUE and associated confidence intervals for an intercept-only model. Those values are included for information but will not be directly used to predict total bycatch. The best model will always be refitted, ensuring the contribution of random effects, if any, is taken into account when predicting on corresponding substrata of the fishing effort.

In order to perform the estimation, the user needs to define the level at which bycatch needs to be assessed (`analysis_resolution` e.g. species, species x region, species x gear, species x gear x region etc), the unit of effort (`effort_term` e.g. days at sea, number of hauls, number of lines etc), the response (`response_term` e.g. number of individuals, number of incidents), as well as any additional variables to use as fixed effect (`fi_terms` e.g. number of hooks, soak time, use of mitigation measures) and/or random effect (`re_terms`, e.g grouping variables, for instance source of bias). Observations can be weighted using custom, user-defined weights (`weights_values`). More detailed descriptions of the argument is available in the report or in the function file itself, accessible through `View(calc_bpue)`.


``` r
#The needle is a table containing all variables necessary to describe the level of assessment. Here, we want one estimate for each combination of region, gear and species. 

needle_a <- unique(bycatch_data[, .(region, gear, species)])

#We can then call the main function, specifying all the required arguments.
bpue1 <- calc_bpue(needle = needle_a, min_re_obs = 2, dat = bycatch_data, response = "n_ind", effort_term = "days", analysis_resolution = c("species", "region", "gear"), re_terms = c("area", "vessel_size", "season"))

#Let's have a look at the results

head(bpue1)
```

```
##    species region     gear observed_bycatch observed_effort
##     <fctr> <fctr>   <fctr>            <int>           <num>
## 1:       A  North Gillnets             2252        3520.078
## 2:       B  North Gillnets              140        3979.041
## 3:       C  North Gillnets               11        3577.826
## 4:       D  North Gillnets              727        3765.654
## 5:       E  North Gillnets              480        3823.950
## 6:       F  North Gillnets               88        3657.777
##                                          model        bpue          lwr
##                                         <char>       <num>        <num>
## 1: resp ~ 1 + (1 | vessel_size) + (1 | season) 0.286684371 0.0273096559
## 2: resp ~ 1 + (1 | vessel_size) + (1 | season) 0.009640416 0.0006907228
## 3:                                    resp ~ 1 0.003314192 0.0015724962
## 4:                                    resp ~ 1 0.176421287 0.1028579963
## 5:                                    resp ~ 1 0.150788339 0.0852343166
## 6:                                    resp ~ 1 0.016972921 0.0055932936
##            upr replicates base_model_heterogeneity alternative_models_flag
##          <num>      <int>                   <lgcl>                  <lgcl>
## 1: 3.009482389         30                     TRUE                    TRUE
## 2: 0.134551260         30                     TRUE                    TRUE
## 3: 0.006984988         30                    FALSE                    TRUE
## 4: 0.302596510         30                     TRUE                   FALSE
## 5: 0.266760197         30                     TRUE                   FALSE
## 6: 0.051504544         30                     TRUE                    TRUE
```

``` r
table(is.na(bpue1$bpue))
```

```
## 
## FALSE  TRUE 
##    45    15
```

We can see that bycatch rate was successfully estimated for 45 cases out of 60. That corresponds to any combination in which at least one individual was reported bycaught:


``` r
table(bpue1$observed_bycatch>0)
```

```
## 
## FALSE  TRUE 
##    15    45
```

What of the 15 remaining cases? We know no bycatch was detected, but can we confidently say there is no bycatch, or is monitoring too low to support such a claim ?

## Accepting BPUE = 0 based on simulations and monitoring effort

The `zero_assessments()` function implements the process used in BEAM to “retrieve” some BPUE values of 0. A zero BPUE is accepted for a given context (e.g. a given species, region and gear combination) when the probability of observing no bycatch if it actually occurred at a “true” rate of X, knowing the observed effort, falls below p = 0.01 (meaning it would be very unlikely that bycatch would have gone undetected).

The function runs on the BPUE table output from `calc_bpue`. It uses as inputs the best model column, to re-examine only cases where no model was built because there were no bycatch records in the data, and the associated observed monitoring effort. The user needs to specify again the resolution at which the analysis was performed (`analysis_resolution`). By default, X, the “true” bycatch probability, is set to 0.001 (1 in 1000 fishing operations), defined as “very rare bycatch” by WKPETSAMP (ICES, 2024) and retained for use in BEAM’s zero assessment during WKBEAM (ICES, 2026).

For cases that meet the \<0.01 probability criterion, the NA *bpue* value is replaced with 0. The lower bound of the confidence interval is also set to 0, while the upper bound is set parametrically as ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAIgAAAAsCAMAAACJ6QCYAAAAAXNSR0IArs4c6QAAANhQTFRFAAAAFT1jFT2AFT2bFWKbFWKAFWK1FYS1FYTORT1jRT2bRT2ARWJjRWKARWKbRWK1RYS1RYTORaXORaXnbz1jbz2bbz2Ab2Jjb2KAb2Kbb2K1b4TOb6Wbb6XOb6Xnb8P/b8PnlWJjlWKAlWKblYSAlYSblYS1laWAlaXnlcO1lcPOlcP/lcPnleHOleH/leHnuYRjuYSAuaWAucObucPnuf//ueH/ueHn3aWA3aWb3cOb3cO13cPn3eG13eHO3f/O3f//3eHn3f/n/8Ob/+G1///O/+HO///npf6bHAAAAAF0Uk5TAEDm2GYAAAAJcEhZcwAADsQAAA7EAZUrDhsAAAAZdEVYdFNvZnR3YXJlAE1pY3Jvc29mdCBPZmZpY2V/7TVxAAADLElEQVRYR+1Xi3bSQBBNAtYELS0FBKuWpsV3xaANhuAjNQn7/3/kndkQHhUylMfxHJnDod2wO3v37t27E8M4xIGBAwMHBg4M/M8MBObuQ8Sv1xF1232nfm/3c0hmUFf/CJC06UvwUp/EMWv0N6jTd+iYp5F0qKBfei7NBu7SBgmKgSQn0Zhb2wokFKZKKr7hEQbNiGGoazGZgikSZlsSBDmwsT8Ii4Q13CYhRpwtrxgKMcJkZIwMPxaPWaPHhOfiIRCIcnONEB/xOpSM29nxVF3TLA94b/uOnc8byJPFs6eGt0g+1hg6Ja0o5dqR8o6gzMQ5m9Hnnox1eBlnQJInYIa+9BHMY2/GOgeEQATESh7bPYKrhDYBkjZqkTFyoLfTtqkdkmINYy2W88oeEyAGpPbopdNJG1BImMuMjHU/ZUAOhOCiwRJR7sQ95Ma6ISE091QQrg2tAsMMELGxbgyE3Dhx6KSMu6RTD2YSskUzR2Jj3QxIfAwBVHsMxLNe0XlRfRhbXoPIjXUzIIWjtbGOX7PlbjvCC2RMnpqWIDEbK9ckMbjTVwCH6lrmY2Jw3Zgtm+hqj23DE+RgY2VaUPeoL1OzC6DrF4IE97rMlU2cPnlXnIeNVbsr3fLpsx6oMC3gCexldCw98Xwg58smQjW6SVYCGT0HWjZW7a6gkP5T7lkENqDpqaLv5VlSTQWdxbIJSenRKiDDDw2amdShgZBacA/RevWNHFZ8XNSu9b4FgpKWWecWnUDkLvn0hG5yeMK4hTsDbbasubKJVlcYMcyEl8ZbQ7sS4nKOj6Kw8u02wG8V/2dYHnj2r6aftntp8xNa2HR8gJSfcAeo6/Oopq5udB2ugWRlk8gdlIv16/cDchkyHByatGVdumW/bx5XsXVBBxpGJ6wM06J17etP9oSHgqDqAJ2yZFzjZ2WTrMgAydpYl75S0BvEeRTX7wAkwKKB6e4HPie/b/kJHTV257AX27g3vi9sg7Qy90pvtcOP3/zd0MACJgwtFA94m0IrcaCUetooQ9T0fqXF3rIuwEH5q7VQOrKhSSJtrFN1SjI+tE+QX34PzbClcft6A/8D7kB2SZeOxoAAAAAASUVORK5CYII=). The function returns an updated version of the BPUE table, which includes the zero values where statistically supported.


``` r
#This will return the same object as the first function (bpue table) but some rows, with sufficient monitoring coverage to statistically support a bycatch rate of 0, will now have 0 in the "bpue" columns. Others remain as NA is no bycatch was reported, or the intitial bpue value if it was computed before.

bpue1_with_zeros <- zero_assessments(bpue = bpue1, analysis_resolution = c("species", "region", "gear") )

#Let's have a look at the results now:
table(is.na(bpue1_with_zeros$bpue))
```

```
## 
## FALSE  TRUE 
##    51     9
```

We can see that a bycatch rate of 0 ind/DaS can be statistically supported in 4 cases, considering the amount of monitoring effort available. If bycatch was occurring for those species, gear and region combinations, even rarely (1 individual in 1000 days at sea), it would have most likely be reported at least once.

# Total Bycatch

## Basic estimation of total bycatch

The `calc_total()` function runs on the BPUE table (`bpue1`) generated in `calc_bpue()`*,* or its updated version, `bpue1_with_zeros`, obtained from `zero_assessments()`. The user also needs to supply the initial monitoring dataset, an additional dataset containing the total fishing effort, and specify again the names of the effort and response columns (`effort_term` and `response_term`), the resolution at which the analysis was performed (`analysis_resolution`), and any weights to be used (`weight_values`) if required. A **filter** can be included in the function call to predict bycatch on a **subset** of the total fishing effort in the table (for instance, only for one métier, or one year, or one season...).


``` r
bycatch1 <- calc_total(bpue = bpue1_with_zeros, dat = bycatch_data, response = "n_ind", effort_term = "days", analysis_resolution = c("species", "region", "gear"), fishing = fishing_data)

#Filter example: Estimating bycatch in winter only (The fishing effort for winter will be used to make predictions, and the predicted numbers will be added to the observed bycatch for winter only)
bycatch_winter <- calc_total(bpue = bpue1_with_zeros, dat = bycatch_data, response = "n_ind", effort_term = "days", analysis_resolution = c("species", "region", "gear"),fishing = fishing_data, filter = list(season = "winter"))
```

Could we estimate total bycatch for all species x gear x region combinations for which we had a BPUE estimate ?


``` r
table(bycatch1$message)
```

```
## 
## More levels available in fishing effort than in monitoring data for at least one random effect - Cannot estimate total bycatch for the levels in which BPUE is unknown 
##                                                                                                                                                                      4 
##                                                                                                                                                                     OK 
##                                                                                                                                                                     46 
##                                                                                                          Unable to explain heterogeneity with candidate random effects 
##                                                                                                                                                                     10
```

Looking at the messages, we can see that, among the 49 BPUE estimates available, total bycatch could not always be estimated:

-   In four cases, at least one random effect influencing BPUE had more levels represented in the total fishing effort data than in the monitoring data. This indicates that part of the fleet was not monitored. For these cases, it is still possible to estimate partial bycatch, considering only the portion of the fleet for which bycatch-rate information is available (see the following section).

-   In another 10 cases, heterogeneity was detected across substrata, but none of the candidate random effects provided sufficient support to improve the model significantly. In these situations, the `alternative_models_flag` indicates whether at least one alternative model containing random effects falls within 2 AIC units of the selected model. Users may then choose to manually refit one of these alternative models in order to estimate total bycatch.

## Cases in which heterogeneity occurs within the fleet, but not all segments are monitored

### Option #1: Estimation of partial bycatch - excluding unmonitored levels

The `calc_partial()` function is similar to the previous one, but runs on `bycatch1`, the output of `calc_total()`, **after** estimation of total bycatch, rather than on `bpue1`. It considers all cases in which total bycatch could not be estimated (due to unmonitored levels for certain variables affecting BPUE), refits the model and predicts bycatch on a subset of the fishing effort, excluding those missing levels. While informative, the resulting figures are **only partial** and are likely an underestimation of the actual amount of bycatch. They should therefore not be considered reliable estimates of bycatch for a given fishery, as not all segments of the fleet are represented. The extent to which the missing segments differ in terms of bycatch, whether similar, higher, or lower, will vary depending on context. A new column in the output table indicates whether the estimate is for total bycatch (unchanged output from previous object) or partial only (new added values).


``` r
bycatch1 <-  calc_partial(tot = bycatch1, dat = bycatch_data, response = "n_ind", effort_term = "days", analysis_resolution = c("species", "region", "gear"),fishing = fishing_data)
```

We now obtain estimates for the four cases in which monitoring coverage was incomplete. Looking at the input data, we can notice that this is because, in all four cases (all involving gillnets), BPUE varied seasonally, while no monitoring data were available for that gear in fall.

### Option #2: Estimation of total bycatch - using the model's average pattern to predict bycatch in unmonitored segment of the fleet

This option is still in development, but will be available shortly. It will allow users to predict bycatch even for levels for which no monitoring is available and therefore BPUE is unknown.

In these cases, predictions will be based on the model's average (population-level) random effects, rather than level-specific estimates. The assumption will be that the unmonitored portion of the fleet behaves similarly to the fleet as a whole and therefore follows the average pattern captured by the model. If the user has pre-existing knowledge of the fishery and suspects that the unmonitored levels are widely different from the rest, partial estimates may remain a better option.

# Reliability of Bycatch Estimates

The `reliability_estimation` function performs and combines two reliability measures:

-   A measure of BPUE sensitivity to exclusion of individual observations, based on an RMSE comparison with leave-one-out jackknife estimates.

-   The breadth of the confidence interval around the total bycatch estimates (estimates are considered uncertain if the CI spans more than 2 orders of magnitude)

It returns a data table identical to the input (`bycatch1`) with new columns for each of the individual and combined criteria. TRUE/FALSE values indicate whether the check was passed or failed. Only estimates passing both checks are considered more reliable.

The function runs on `bycatch1`. It does not assess whether the estimate is total or partial. but because partial estimates should generally be regarded as uncertain, the reliability output can be manually overwritten for those specific cases.


``` r
#Hard code that if the partial option has been used, they should be uncertain. 
# Add the monitoring coverage reliability test ? Threshold of <0.001. Maybe as an argument

bycatch1_reliability <- reliability_estimation(bycatch1, dat = bycatch_data, effort_term = "days", response_term = "n_ind", analysis_resolution = c("species", "region", "gear"))

#If partial estimates have been generated using calc_partial, the estimate_quality column will exist in the table, we can use it to overwrite the overall reliability value

bycatch1_reliability[estimate_quality == "Partial estimate only", Overall_reliability := FALSE]
```

# Visualisation

The code below provides a few example of visualisation for the results of this specific dataset.

## Success / Failure matrix plots

Matrix plots can be used to summarise the overall outcome of BEAM across species, region and gear. This can give a quick idea of how much is known or remains unknown in a given fishery. Here, cells are colored differently whether total bycatch was estimates, partial, or none. The outline of the cell indicates whether that estimate can be considered reliable or whether it's uncertain.


``` r
#This is an example of data visualisation, showing estimate success rate, reliability, overview of the fishery.

bycatch1_reliability[is.na(tot_mean), estimate_quality := NA]
bycatch1_reliability[tot_mean == 0, Overall_reliability := TRUE]


ggplot(bycatch1_reliability, aes(x = paste("Species",species), fill = as.factor(estimate_quality),col = as.factor(Overall_reliability))) +  geom_bar(position = "fill", width = .7, lwd=1) + coord_flip() + labs(y = "", x = "", fill = "", col = "")  +  scale_color_manual(
    values = c("#DD6E41","#129484"),
    drop = FALSE,
    na.value = "grey95",
    labels = c("Uncertain estimate", "More reliable estimate")
  )+ scale_fill_manual(
  values = c(
    "Total estimate" = "#129484",
    "Partial estimate only" = "#796E62"),
    na.value = "grey95",  drop = FALSE)+ theme(strip.background = element_rect(fill = "white"),panel.background = element_rect(fill = "grey95"),axis.title = element_text(size = 10), axis.text = element_text(size = 7)) + facet_grid(str_to_title(region) ~ toupper(gear), drop = FALSE) + scale_x_discrete(labels = str_to_title)+ scale_y_continuous(labels = NULL, breaks = NULL)+theme_minimal()+theme(legend.position = "right", legend.box="vertical")+ guides(
               colour = guide_legend(
                 override.aes = list(fill = "white")))
```

![](CIBBRiNA_bycatch_estimation_toolkit_implementation_example_github_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

## BPUE plots

The plot below shows the bycatch rate (BPUE) for each species in each gear and region, when it could be estimated. The first panel includes bycatch rates estimated from the data, and the second panel the BPUE accepted to be 0, given the monitoring effort.


``` r
#Plots of BPUEs

bpue1_with_zeros$gear <- toupper(bpue1_with_zeros$gear)
bpue1_with_zeros$gear <- as.factor(bpue1_with_zeros$gear) #Make sure gear is a factor for facets to work properly.

plot1 <- ggplot(data=bpue1_with_zeros[bpue>0,], aes(x= paste("Species",str_to_title(species)), y= bpue, col = str_to_title(region))) + geom_point(position=position_dodge(.5), size=2, alpha=.8)+ geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2,         position=position_dodge(.5)) +scale_color_manual(values= c("#129484", "#DD6E41"))+scale_y_log10(labels = label_log())+coord_flip()+theme_minimal()+facet_grid(~gear, scales = "free_y", space = "free_y")+labs(y="BPUE (n_ind/DaS) - Log10", x="", title="Modelled BPUEs", col="Region")+theme(panel.background = element_rect(colour = "#796E62", fill =NA ), legend.position = "none")

plot2 <- ggplot(data=bpue1_with_zeros[bpue==0,], aes(x= paste("Species",str_to_title(species)), y= bpue, col = str_to_title(region))) + geom_point(position=position_dodge(.5), size=2, alpha=.8)+ geom_errorbar(aes(ymin=lwr, ymax=upr), width=.2,         position=position_dodge(.5)) +scale_color_manual(values= c("#129484", "#DD6E41"))+coord_flip()+theme_minimal()+facet_grid(~gear, scales = "free_y", space = "free_y", drop=FALSE)+labs(y="BPUE (n_ind/DaS)", x="", title="BPUEs = 0 - Accepted from simulations", col="Region")+theme(panel.background = element_rect(colour = "#796E62", fill =NA ), legend.position = "bottom")+scale_y_continuous(labels = scales::label_number(accuracy = 0.0001))


plots_BPUE <- ggpubr::ggarrange(plot1, plot2, nrow = 2, heights = c(2,1.3))


plots_BPUE
```

![](CIBBRiNA_bycatch_estimation_toolkit_implementation_example_github_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

## Total bycatch estimate plots, with reliability information

The following plot presents estimates of total bycatch, where available, for each species, gear and region combinations, along with information on the reliability of the estimates.


``` r
bycatch1_reliability$gear <- toupper(bycatch1_reliability$gear)
bycatch1_reliability$gear <- as.factor(bycatch1_reliability$gear) ##Make sure gear is a factor for facets to work properly.

plot3 <- ggplot(data=bycatch1_reliability[tot_mean>0 & estimate_quality == "Total estimate",], aes(x= paste("Species",str_to_title(species)), y= tot_mean, col = str_to_title(region), pch = Overall_reliability)) + geom_point(position=position_dodge(.5), size=2, alpha=.8)+ geom_errorbar(aes(ymin=tot_lwr, ymax=tot_upr), width=.2, position=position_dodge(.5)) +scale_color_manual(values= c("#129484", "#DD6E41"))+scale_y_log10(labels = label_log())+coord_flip()+theme_minimal()+facet_grid(~gear, scales = "free_y", space = "free_y")+labs(title = "Predicted Total Bycatch", y="Total bycatch (n_ind) - Log10", x="", col="Region", shape = "Reliability")+theme(panel.background = element_rect(colour = "#796E62", fill =NA )) +scale_shape_manual(values = c(4,16),labels = c("Uncertain", "More reliable"))+theme(legend.position = "top")

plot4 <- ggplot(data=bycatch1_reliability[tot_mean>0 & estimate_quality == "Partial estimate only",], aes(x= paste("Species",str_to_title(species)), y= tot_mean, col = region)) + geom_point(position=position_dodge(.5), size=2, alpha=.8, pch=4)+ geom_errorbar(aes(ymin=tot_lwr, ymax=tot_upr), width=.2, position=position_dodge(.5)) +scale_color_manual(values= c("#129484", "#DD6E41"))+scale_y_log10(labels = label_log())+coord_flip()+theme_minimal()+facet_grid(~gear, scales = "free_y",  drop = FALSE)+labs(title = "Predicted Partial Bycatch",y="Partial bycatch (n_ind) - Log10", x="", col="Region")+theme(panel.background = element_rect(colour = "#796E62", fill =NA ))+theme(legend.position = "none")

plot5 <- ggplot(data=bycatch1_reliability[tot_mean==0,], aes(x= paste("Species",str_to_title(species)), y= tot_mean, col = region)) + geom_point(position=position_dodge(.5), size=2, alpha=.8)+scale_color_manual(values= c("#129484", "#DD6E41"))+coord_flip()+theme_minimal()+facet_grid(~gear, scales = "free_y",  drop = FALSE)+labs(y="Total bycatch (n_ind)", x="", col="Region")+theme(panel.background = element_rect(colour = "#796E62", fill =NA ))+theme(legend.position = "none")+scale_y_continuous(breaks = c(0,1))


bottom_plots <- ggpubr::ggarrange(plot5, plot4, ncol = 1,heights = c(1, 1.2))

plots_bycatch <- ggpubr::ggarrange(plot3, bottom_plots, ncol = 1, heights = c(1.3, 1))


plots_bycatch
```

![](CIBBRiNA_bycatch_estimation_toolkit_implementation_example_github_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

![](C:/Users/mpommier/OneDrive%20-%20Marine%20Institute/Morgane/CIBBRiNA/Communication%20materials/Funding%20banner.png)
