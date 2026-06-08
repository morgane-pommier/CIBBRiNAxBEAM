
#' Calculates total bycatch for bycatch models fitted in BEAM
#' 
#' @description
#' `calc_total` calculates total bycatch based on a BPUE.
#' 
#' #' @details
#' If there are multiple rows in needle, each row is processed separately, and
#' combined into one data.table with nrow equal to the number of rows in needle.
#' 
#' @param bpue data.table output from the BPUE estimation function (calc_bpue.r), with values for for `analysis_resolution`, used to subset observations from `dat`. 
#' @param response Character string. Name of the response variable. Typically, a column containing the number of individuals or bycatch incidents (positive integer values).
#' @param effort_term Character string. Name of the variable to be used as an offset (e.g. number of days at sea, number of hooks, number of nets …), in both `dat` and `fishing`. By default, the log of this covariate will be used as offset in the model. 
#' @param analysis_resolution vector of column names, whose unique combinations are used to split the data given in `dat` and `fishing`. Must be the same vector used in calc_bpue.r. Note that if `species` is among those columns, it is not used when subsetting `fishing` (because the fishing effort dataset does not include information about species).
#' @param dat data.table with observed (monitored) fishing effort data. Same `dat` input table used in calc_bpue.r 
#' @param fishing data.table with all fishing effort data.
#' @param verbose Logical, passed on to predict_response. Set to FALSE if nrow(bpue) > 1.
#' @param filter optional named list for subsetting `fishing` and the observed bycatch summation (but NOT model fitting. BPUE is estimated based on all data available, weighted is applicable). Each element name must be a column in both `fishing` and `dat`, and each element value is a vector of levels to retain. E.g. `list(year = 2023, ecoregion = "North Sea")`. Default NULL
#' @param include.weights boolean indicator indicating whether observations from the five most recent years should have double weight as compared to older data. 
#' @param weight_values Optional. Numeric vector of observation weights. Must be of the same length that dat. Only used if include.weights = TRUE. Default is NULL.
#' @returns A data.table with all columns given in `analysis_resolution`, and additional columns showing the model formula, total bycatch estimate, lower and upper confidence intervals, and total fishing effort. See details.
#' @seealso [calc_bpue()]
#' @export

calc_partial <- function(tot, analysis_resolution, dat, fishing, verbose = TRUE, include.weights=FALSE,response, effort_term,weights_values, filter=NULL) {
  
  # parallelization support
  if (nrow(tot) > 1) {
    ret <- foreach(i = 1:nrow(tot), 
                   .export = "calc_total", # <- not 100% sure this line is needed.
                   .final = rbindlist,
                   .packages = c("data.table", "glmmTMB", "emmeans", "ggeffects")) %dopar% {
                     calc_total(tot = tot[i], analysis_resolution = analysis_resolution, dat = dat, fishing = fishing, verbose = FALSE, include.weights=FALSE,response, effort_term = effort_term,weights_values, filter)
                   }
    return(ret)
  }
  
  #Apply any filter to fishing if needed. Checks that the filtering columns are present in both fishing effoert and monitoring data, so the observed bycatch can be filtered accordingly before summation at the end
  
  filter_fishing <- NULL
  filter_dat     <- NULL
  
  if (!is.null(filter)) {
    if (!is.list(filter) || is.null(names(filter))) {
      stop("`filter` must be a named list, e.g. list(year = 2023, ecoregion = 'North Sea')")
    }
    missing_in_fishing <- names(filter)[!names(filter) %in% names(fishing)]
    missing_in_dat     <- names(filter)[!names(filter) %in% names(dat)]
    if (length(missing_in_fishing) > 0)
      warning(paste("filter column(s) not found in `fishing` and will be ignored:", 
                    paste(missing_in_fishing, collapse = ", ")))
    if (length(missing_in_dat) > 0)
      warning(paste("filter column(s) not found in `dat` and will be ignored for observed bycatch summation:", 
                    paste(missing_in_dat, collapse = ", ")))
    # only apply valid columns
    filter_fishing <- filter[names(filter) %in% names(fishing)]
    filter_dat     <- filter[names(filter) %in% names(dat)]
  }
  
  #Filter fishing data if needed
  fishing_filtered <- fishing
  if (!is.null(filter) && length(filter_fishing) > 0) {
    for (col in names(filter_fishing)) {
      fishing_filtered <- fishing_filtered[get(col) %in% filter_fishing[[col]]]
    }
  }
  
  #Focus only on cases which were not calculated due to missing levels in monitoring
  
  tot <- tot[message == "More levels available in fishing effort than in monitoring data for at least one random effect - Cannot estimate total bycatch for the levels in which BPUE is unknown"],
  
  #Create standardised column to use in data manipulation and modelling step further below. Avoids having to use get() of dat[[]] throughout the entire code.
  dat[, resp   := get(response)]
  dat[, effort := get(effort_term)]
  if (isTRUE(include.weights)) dat[, weights := get(weights_values)]
  
  dat <- dat[tot, on = analysis_resolution]
  dat <- dat[effort > 0]
  
  dat[, (analysis_resolution) := lapply(.SD, as.factor), .SDcols = analysis_resolution]
  dat[, logEffort := log(effort)]
  dat[, resp := as.integer(resp)] 
  
  # add observation weights to data(if monitoring within the last five years 
  # use full weight, if older data weight observations half as strongly in the likelihood)
  #dat[, weights := ifelse(year >= (max(years)-4),1,0.5)] OLD VERSION, KEEP AS AN EXAMPLE
  
  ret <- copy(tot)
  
  form <- as.formula(tot$model)
  re <- lme4::findbars(form) # random effects part of model formulation (if any)
  re <- sapply(re, function(x) as.character(x[[3]]))
  re.n <- length(re)
  
  
  # sum up total fishing effort per combination of all levels of the random effects
  fishing_filtered <- fishing[tot[, ..analysis_resolution],
                              on = analysis_resolution[analysis_resolution != "species"],
                              nomatch = 0]
  fishing_filtered[, effort := get(effort_term)]
  
  # for each random effect, keep only levels present in monitoring data
  for (r in re) {
    if (r %in% names(fishing_filtered) && r %in% names(dat)) {
      monitored_levels <- unique(dat[[r]])
      fishing_filtered <- fishing_filtered[get(r) %in% monitored_levels]
    }
  }
  
  if (nrow(fishing_filtered) == 0) {
    ret$message <- "No fishing effort remaining after restricting to monitored levels"
    return(ret)
  }
  
  #Create the partial object
  tot_partial <- fishing_filtered[, .(effort = sum(effort)), by = re]
  tot_partial <- tot_partial[complete.cases(tot_partial)]

  # Refit best model
  if(isTRUE(include.weights)){
    best <- glmmTMB(formula = form, offset = logEffort, family = nbinom2, data = dat, weights = weights)
  } else {
    best <- glmmTMB(formula = form, offset = logEffort, family = nbinom2, data = dat)  
  }
  
  dat_filtered <- dat
  if (!is.null(filter) && length(filter_dat) > 0) {
    for (col in names(filter_dat)) {
      dat_filtered <- dat_filtered[get(col) %in% filter_dat[[col]]]
    }
  }
  
    
    # sum monitoring effort, and number of bycaught animals, per random factor 
    tot_obs <- dat_filtered[,(re) := lapply(.SD, as.factor), .SDcols = re]
    tot_obs <- tot_obs[,
                       .(observed_effort = sum(effort), 
                         observed_bycatch = sum(resp)), 
                       by = re]
    
    tot_partial[, (re) := lapply(.SD, as.factor), .SDcols = re] # Do we treat YEAR correctly here? or does year need to be exempted?
    
    # join total fishing effort and monitored effort 
    tot_partial <- tot_obs[tot_partial, on=re]
    tot_partial[is.na(observed_effort), observed_effort := 0] # replace na with 0
    tot_partial[is.na(observed_bycatch), observed_bycatch := 0] # replace na with 0
    tot_partial[,unmonitored_fishing_effort := effort - observed_effort] #
    tot_partial[, logEffort := log(unmonitored_fishing_effort)]
    
    
    pred <- lapply(1:nrow(tot), function(i) {
      type_arg <- ifelse(packageVersion("ggeffects") >= "2.0.0", "random", "re")
      
      p <- ggpredict(model = best,
                     terms = tot_partial[i, ..re],
                     condition = c(logEffort = tot_partial$logEffort[i]),
                     type = type_arg,
                     interval = "confidence",
                     verbose = verbose)
      
      p <- as.data.frame(p) 
      data.table(mean = p$predicted,
                 lwr = ifelse(!is.null(p[["conf.low"]]), p[["conf.low"]], NA_real_),
                 upr = ifelse(!is.null(p[["conf.high"]]), p[["conf.high"]], NA_real_))
    }) |> rbindlist()
    
    ret[, c("tot_mean", "tot_lwr", "tot_upr") := as.list(colSums(pred) + sum(tot_partial$observed_bycatch))] # prediction for unmonitored fishing effort + observed bycatch in monitoring
    ret$fishing_effort <- sum(tot_partial$effort) # retain total fishing effort or total unmonitored fishing effort?
    
  
  return(ret)
}

