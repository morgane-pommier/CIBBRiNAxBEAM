## ---------------------------
## Authors :  WGBYC ToR C members. Original scripts available at https://github.com/dlusseau/BEAM/tree/main
## Adapted for CIBBRiNA users by MP.
## ---------------------------
#' 
#' @description
#' `reliability_estimation` Runs on the output of `calc_total()`. Fits a full glmmTMB model and leave-one-out jackknife models for each unique bycatch estimate. Computes RMSE between BPUE predictions from both models on the log scale. Returns the input data (total estimates) enriched with uncertainty and reliability metrics combining 1) model stability across jackknife folds at the BPUE level with 2) the breadth of the confidence interval around the total bycatch estimate. The bycatch estimate can be considered reliable is both checks are passed (overall_reliability).
#' 
#' @param tot data.table output from the total bycatch estimation function `calc_total()`, with values for for `analysis_resolution`, used to subset observations from `dat`. Also contains a column for `model` and `tot_mean`, the bycatch estimate value, if any.
#' #' @param dat data.table with observed (monitored) fishing effort data. Same `dat` input table used in `calc_bpue()` and `calc_total()`
#' @param response Character string. Name of the response variable. Typically, a column containing the number of individuals or bycatch incidents (positive integer values).
#' @param effort_term Character string. Name of the variable to be used as an offset (e.g. number of days at sea, number of hooks, number of nets …), in `dat`. By default, the log of this covariate will be used as offset in the model. 
#' @param analysis_resolution vector of column names, whose unique combinations are used to split the data given in `dat`. Must be the same vector used in `calc_bpue()` and `calc_total()` 
#' @param include.weights boolean indicator indicating whether observations from the five most recent years should have double weight as compared to older data. Must match the value used in `calc_bpue`. 
#' @param weight_values Optional. Numeric vector of observation weights. Must be of the same length that dat. Only used if include.weights = TRUE. Default is NULL.
#' @returns A data.table identical to tot, with additional columns: 
#' - rmse-val: root mean square error of jackknife vs full model predictions (log scale)
#' - factor:  RMSE back-transformed from log scale 
#' - lower_factor_prop — lower bound of the 68% uncertainty interval as a percentage: 100(1 - 1/factor)
#' - upper_factor_prop — upper bound of the 68% uncertainty interval as a percentage: 100(factor - 1)
#' - reliability_rmse —  TRUE/FASLE, whether both factor proportions are ≤ 25%
#' - CI_breadth_<2 — TRUE/FASLE, whether the breadth of the confidence interval for the total bycatch estimate is smaller than 2 orders of magnitude.
#' - Overall_reliability — TRUE/FASLE, whether both checks above are passed, or at least one is failed.
#' @seealso [calc_bpue(), calc_total()]
#' @export

#Difference with BEAM original code: BEAM reliability_estimation_p.r scripts reads data, calculates RMSE using a calc_rmse function, and then fills in new columns in the input data. Here, everything is happening in a single reliability_estimation.r function, and input data is supplied in the function call. The updated function also includes the breath of CI calculations and assessment.

reliability_estimation <- function(tot,
                      dat,
                      analysis_resolution,
                      effort_term,
                      include.weights = FALSE,
                      weights_values  = NULL) {
  

#Parallel set-up, updated with CIBBRiNA functions arguments
  
  if (nrow(tot) > 1) {
    
    if (exists("BEAM_pb")) {
      BEAM_pb$terminate()
      rm("BEAM_pb", envir = .GlobalEnv)
    }
    BEAM_pb <<- progress_bar$new(
      format = "reliability_estimation :percent :current/:total [:bar] :elapsed | eta: :eta",
      total = nrow(tot),
      width = 60)
    
    opts <- list(progress = BEAM_progress)
    
    ret <- foreach(i = 1:nrow(tot),
                   .export = "reliability_estimation",
                   .final = rbindlist,
                   .packages = c("data.table", "glmmTMB", "ggeffects", "emmeans"),
                   .options.snow = opts) %dopar% {
                     reliability_estimation(
                       tot = tot[i],
                       dat = dat,
                       analysis_resolution = analysis_resolution,
                       effort_term = effort_term,
                       include.weights = include.weights,
                       weights_values = weights_values
                     )
                   }
    BEAM_pb$terminate()
    rm("BEAM_pb", envir = .GlobalEnv)
    return(ret)
  }
  
  
  #Create new columns for RMSE calculations:
  
  ret <- copy(tot)
  ret[, c("RMSE", "factor", "lower_factor_prop", "upper_factor_prop",
          "reliability_rmse", "CI_breadth_pass", "Overall_reliability",
          "tot_lwr_log", "tot_upr_log", "CI_breadth") :=
        list(NA_real_, NA_real_, NA_real_, NA_real_, NA, NA, NA,
             NA_real_, NA_real_, NA_real_)]
  
  #Skip cases that do not have a total bycatch estimate or no usable model. Partial estimates will by default have a unreliable estimate. Do we want to change that and still compute RMSE on partial models ?
  
  if (is.na(tot$model) || tot$model %in% c("none", "only one") || is.na(tot$tot_mean)) {
    return(ret)
  }
  
  
    # helper function (to avoid repeated code). Same as original BEAM, just added weights option
  fit_and_predict <- function(form, data, levels = NULL) {
    fit <- suppressMessages(tryCatch({
      if (isTRUE(include.weights)) {
        glmmTMB(formula = form, offset = logEffort, family = nbinom2, data = data,weights = weights)
      } else {
        glmmTMB(formula = form, offset = logEffort,family = nbinom2,data = data)
      }
}, error = function(e) NULL))
    re <- lme4::findbars(form) # random effects part of model formulation (if any)
    re <- sapply(re, function(x) as.character(x[[3]]))
    re.n <- length(re) # number of random effects
    
    if (re.n == 0) {
      pred <- as.data.table(emmeans(fit, ~1, type = "response", offset = log(1)))
      pred[, predicted := response]
         } else {
             type_arg <- ifelse(packageVersion("ggeffects") >= "2.0.0", "random", "re")
             pred <- ggpredict(model = fit, type =  type_arg, terms = re, verbose = FALSE)
             pred <- as.data.table(pred)
             
             if (!"facet" %in% colnames(pred)) {
                 pred[, facet := "dummy"]
             }
             pred <- pred[order(x, group, facet)]
         }
     
         list(fit = fit, re = re, re.n = re.n, pred = pred)
    }
    
    needle <- tot[, ..analysis_resolution]
    data_subset <- dat[needle, on = analysis_resolution, nomatch =0]
    
    #Get column names for later steps, and compute logEffort
    
    data_subset[, effort := get(effort_term)]
    if (isTRUE(include.weights)) data_subset[, weights := get(weights_values)]
    data_subset <- data_subset[effort > 0]
    data_subset[, (analysis_resolution) := lapply(.SD, as.factor), .SDcols = analysis_resolution]
    data_subset[, logEffort := log(effort)]
    
    # need at least 2 rows for LOO to be meaningful
    if (nrow(data_subset) < 2) return(ret)

    # Extract the formula of the model
    form <- as.formula(tot$model)
    if (is.null(form)) return(ret)
    
    # Fit full model
    full_model <- fit_and_predict(form, data_subset)
    if (is.null(full_model)) return(ret)

    #LOO jackknife 
    
    diffs <- sapply(1:nrow(data_subset), function(j) {
        jack <- fit_and_predict(form, data = data_subset[-j,])
        log_diff <- log1p(full_model$pred$predicted) - log1p(jack$pred$predicted) #What happens if there are less predictions in the jackknife model ? Because a level of a random effect might have been dropped ? Is that a possibility ?
        n <- ifelse(jack$re.n == 0, 1, full_model$re.n) #Why are we dividing by number of random effects in the full model ? Can it be a different value that the jackknife model ?
        sum(log_diff^2, na.rm = TRUE) / n
    })
    
    #This was happening after calc_rmse in the original function, brought back everything in a single function here.
    rmse_val <- sqrt(sum(diffs, na.rm=TRUE)/(nrow(data_subset)-1)) #Correction: Added -1 here to match description at top of the function.
    # Transform back from log scale
    factor_val <- exp(rmse_val)
    # The lower limit of the 68% interval is: 100(1-1/factor) # in percentage
    lower_factor_prop_val <- 100 * (1 - 1 / factor_val)
    # The upper limit of the 68% interval is: 100(factor−1) # in percentage
    upper_factor_prop_val <- 100 * (factor_val - 1)
    # RELIABILITY LIMIT ON BOTH CONFIDENCE INTERVALS OF FACTOR
    # Take the 25% limit of the lower % CI of factor (=80% similarity)
    # Take the 25% limit of the upper % CIs of factor (=80% similarity)
    # Which are below 25%
    
    ret[, `:=`(
      tot_lwr_log = log10(tot_lwr + 1),
      tot_upr_log = log10(tot_upr + 1),
      CI_breadth  = log10(tot_upr + 1) - log10(tot_lwr + 1)
    )]
    
    ret[, c("RMSE", "factor", "lower_factor_prop", "upper_factor_prop") :=
          list(rmse_val, factor_val, lower_factor_prop_val, upper_factor_prop_val)]
    
    #Make all the TRUE/FALSE checks
    ret[, reliability_rmse        := lower_factor_prop <= 25 & upper_factor_prop <= 25]
    ret[, CI_breadth_pass := CI_breadth < 2]
    ret[, Overall_reliability := reliability_rmse & CI_breadth_pass]

    tot_reliability <- ret
    return(tot_reliability)
}






