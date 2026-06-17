## ---------------------------
## Authors :  WGBYC ToR C members. Original scripts available at https://github.com/dlusseau/BEAM/tree/main
## Adapted for CIBBRiNA users by MP.
## ---------------------------

#' @description
#' `calc_bpue` uses monitored fishing effort data to generate a number of models of varying
#' complexity to estimate bycatch rate. It then compares those models using the Akaike Information Criterion (AIC), combined with a stability 
#' filter and a parsimony rule
#' @details
#' If there are multiple rows in `needle`, each row is processed separately, and
#' rbind'ed into one data.table with nrow equal to the number of rows in `needle`.
#' If a parallel backend is available, the code will execute in parallel, allowing
#' faster computations.

#' The model selection uses **AIC** with a **Parsimony Rule**:
#' 1. **Model Selection**: Valid combinations of random effects are fitted. Models with 
#'    "Singular Fit" (near-zero variance or perfect correlation) are 
#'    discarded to ensure numerical stability and prevent over-parameterization.
#' 2. **Tie-breaker**: Among stable candidates, if the difference in AIC is < 2, 
#'    the simpler model (fewer random effects) is chosen, but flagged under "alternative_models_flag".

#' @param needle data.table with values for `analysis_resolution`, used to subset observations from `dat`. Any additional columns in `needle` are disregarded.
#' @param analysis_resolution columns, whose unique combinations are used to split the data given in `dat`. By default, equal to colnames(needle)
#' @param response Character string. Name of the response variable. Typically, a column containing the number of individuals or bycatch incidents (positive integer values).
#' @param effort_term Character string. Name of the variable to be used as an offset (e.g. number of days at sea, number of hooks, number of nets …). By default, the log of this covariate will be used as offset in the model. 
#' @param min_re_obs Integer specifying the minimum number of levels needed to include a term as a random effect
#' @param dat data.table with monitored fishing effort data (e.g. since 2017)
#' @param re_terms Optional. Character vector of column names to be considered as candidate random effects (e.g. c("year", "vessel", "area")). Should represent grouping factors that introduce variability or stratification in the data (e.g. sampling strata, fishing vessels, spatial or temporal clusters). Terms with fewer than min_re_obs unique levels will be automatically dropped. Set to NULL to fit without random effects. Defaults to NULL.
#' @param fi_terms Optional. Character vector of column names specifying fixed effect terms. Typically represent additional covariates that may influence bycatch rates but do not define grouping structures in the data (e.g. soak time, net length). Set to NULL to fit models with no fixed effect other than the intercept. Default to NULL. Note: automated model selection is performed over random effects structure only. If comparing different fixed effect specifications, run the function separately for each combination.
#' @param years vector with integers indicating years of assessment. 
#' @param include.weights boolean indicator indicating whether observations from the five most recent years should have double weight as compared to older data. 
#' @param weight_values Optional. Numeric vector of observation weights. Must be of the same length that dat. Only used if include.weights = TRUE. Default is NULL.
#' @returns A data.table with one row for each row in `needle` and all columns given in `analysis_resolution`. Additional columns include, for each estimate: the observed number of individual or incidents (response) observed, the monitoring effort (unit = unit of analysis), the number of replicates (i.e. strata) available to fit the model, the output of the heterogeneity test (TRUE/FALSE) the best model formula and if any, the population level BPUE estimate from an intercept-only model, with associated lower and upper confidence intervals. (I^2). See details.
#' @seealso [calc_total()]
#' @export
#' 
BEAM_progress <- function(n) {
    BEAM_pb$tick(tokens = list(step = n))
}

calc_bpue <- function(needle, analysis_resolution = colnames(needle), min_re_obs = 2, dat, response, years, effort_term, re_terms = NULL, fi_terms=NULL, include.weights = FALSE, weights_values = NULL) {

    t_start <- Sys.time()
    
    # parallelization support
    if (nrow(needle) > 1) {
        
        if (exists("BEAM_pb")) {
            BEAM_pb$terminate()
            rm("BEAM_pb", envir = .GlobalEnv)
			 }
        
        BEAM_pb <<- progress::progress_bar$new(
            format = "calc_bpue :percent :current/:total [:bar] :elapsed | eta: :eta",
            total = nrow(needle),
            width = 60)
        
        opts <- list(progress = BEAM_progress)
        
        ret <- foreach(i = 1:nrow(needle),
                       .export = "calc_bpue", # <- not 100% sure this line is needed.
                       .final = rbindlist,
                       .packages = c("data.table", "glmmTMB", "metafor", "emmeans"),
                       .options.snow = opts) %dopar% {
            calc_bpue(needle = needle[i], analysis_resolution = analysis_resolution, min_re_obs = min_re_obs, response, re_terms=re_terms, fi_terms = fi_terms, effort_term = effort_term, dat = dat, include.weights, weights_values)
                       }
        
        BEAM_pb$terminate()
        rm("BEAM_pb", envir = .GlobalEnv)
        return(ret)
    }

    if (!inherits(dat, "data.table")) {
      warning("dat is not a data.table, converting automatically.")
      dat <- as.data.table(dat)
    }
    
    dat <- dat[needle[, ..analysis_resolution], on = analysis_resolution, nomatch = 0]
    
    #Create standardised column to use in data manipulation and modelling step further below. Avoids having to use get() of dat[[]] throughout the entire code.
    dat[, resp   := get(response)]
    dat[, effort := get(effort_term)]
    if (isTRUE(include.weights)) dat[, weights := get(weights_values)] 
    
    dat <- dat[effort > 0]
    dat[, (analysis_resolution) := lapply(.SD, as.factor), .SDcols = analysis_resolution]
    dat[, logEffort := log(effort)]
	  dat[, resp := as.integer(resp)] 
	  
	  
    # add observation weights to data(if monitoring within the last five years 
    # use full weight, if older data weight observations half as strongly in the likelihood)
    
	  #dat[, weights := ifelse(year >= (max(years)-4),1,0.5)] OLD VERSION, KEEP AS AN EXAMPLE
    
    ret <- needle[, ..analysis_resolution]
    ret[, c("observed_bycatch", "observed_effort", "model", "bpue", "lwr", "upr", "replicates", "base_model_heterogeneity", "alternative_models_flag") :=
            list(sum(dat$resp), sum(dat$effort), "none", NA_real_, NA_real_, NA_real_, nrow(dat), NA, FALSE)]

    if ((nrow(dat) == 0) || (sum(dat$resp, na.rm = TRUE) == 0) || any(dat$resp < 0)) {
        return(ret)
    }
    
    # cop out when we only have 1 row of data
    if (nrow(dat) == 1) {
        bpue <- dat$resp / dat$effort
        lwr <- bpue -1.96 * sqrt(dat$resp / dat$effort^2)
        upr <- bpue +1.96 * sqrt(dat$resp / dat$effort^2)
        ret[, c("bpue", "lwr", "upr", "model") := list(bpue, lwr, upr, "only one")]
        return(ret)
    }
    
    # define base model formula either including or excluding fixed effects terms
    if (is.null(fi_terms)) {
      base_model_formula <- resp ~ 1
    } else {
      base_model_formula <- as.formula(paste("resp ~", paste(fi_terms, collapse = " + ")))
    }
    
    # fit model either including or excluding likelihood weights
    if (isTRUE(include.weights)) {
      base_model <- tryCatch(
        glmmTMB(base_model_formula, offset = logEffort, family = nbinom2, data = dat, weights = weights),
        error = function(e) e$message
      )
    } else {
      base_model <- tryCatch(
        glmmTMB(base_model_formula, offset = logEffort, family = nbinom2, data = dat),
        error = function(e) e$message
      )
    }
    
    #Heterogeneity test using metafor. Should we include the fixed effects here ? 
    heterogeneity <- tryCatch((rma.glmm(xi = resp, ti = effort, measure = "IRLN", data = dat)$QEp.Wld<0.05), error = function(e) e$message)
    
    
    if (!inherits(heterogeneity, "character")) {
        ret$base_model_heterogeneity <- heterogeneity
    }
    
    best <- base_model # if we have less than 5 rows of monitoring data, fit a simple model and we're done.

    # if we have more than 5 rows of monitoring data, fit all possible
    # combinations of the base model and one or more of the terms in the
    # vector below, added to the model as random effects.
    if (nrow(dat) >= 5) {

        # but only consider r.e. terms where the number of unique values 
        # (i.e. levels) is greater than min_re_obs. Note that this effectively
        # prevents any variables specified in the analysis_resolution parameter from being included
        # as random effects in any models, since they will always have length=1.
        re_terms <- re_terms[sapply(re_terms, function(x) length(unique(dat[[x]]))) >= min_re_obs]
        re.i <- do.call(CJ, replicate(length(re_terms), c(TRUE, FALSE), simplify = FALSE)) # all possible combinations of random effect levels
        
        # fixed-effects part of the model
        if (is.null(fi_terms)) {
          fixed_effects <- "1" #set to intercept only model if we do not have fixed effects supplied 
        } else {
          fixed_effects <- paste(fi_terms, collapse = " + ")
        }

        re_candidates <- apply(re.i, 1, function(i) {
          if (all(i == FALSE)) return(fixed_effects)  
          re_part <- sprintf("(1|%s)", re_terms[unlist(i)])
          return(paste(fixed_effects, "+", paste(re_part, collapse = " + ")))  
        })						

	## Integration of two versions of calc_bpue:
	# New model selection function written by KMB, including singularity checks and parsimony rule, 
	# New option of using weights in models, developed by TS

	  candidates <- lapply(re_candidates, function(f_str) {
      curr_formula <- as.formula(sprintf("resp ~ %s", f_str))
      fit <- suppressMessages(tryCatch({
		 if(isTRUE(include.weights)){
             m <- glmmTMB(formula = curr_formula, offset = logEffort, family = nbinom2, data = dat, weights = weights)
            }else{
             m <- glmmTMB(formula = curr_formula, offset = logEffort, family = nbinom2, data = dat)
            }
       
        # Check for Singularity in glmmTMB
        vc <- VarCorr(m)$cond
        is_singular <- FALSE
        for (comp in vc) {
          if (any(diag(comp) < 1e-6)) is_singular <- TRUE # Variance near zero
          if (nrow(comp) > 1) {
            cor_mat <- attr(comp, "correlation")
            if (any(abs(cor_mat[lower.tri(cor_mat)]) > 0.99)) is_singular <- TRUE # Near perfect correlation
          }
        }
        
        if (!is_singular) return(m) else return("Singular") 
      }, error = function(e) e$message))
      return(fit)
    })
    
    # Identify Valid glmmTMB objects
    valid_indices <- which(sapply(candidates, function(x) inherits(x, "glmmTMB")))
    
    if (length(valid_indices) > 0) {
      # Use AIC instead of BIC for more balanced complexity selection
      valid_aics <- sapply(candidates[valid_indices], AIC)
      valid_complexity <- lengths(regmatches(re_candidates[valid_indices], gregexpr("\\+", re_candidates[valid_indices])))
      
      sel_table <- data.table(
        Index = valid_indices,
        AIC = valid_aics,
        Complexity = valid_complexity
      )
      
      sel_table <- sel_table[is.finite(AIC)]
      
      if (nrow(sel_table) > 0) {
        min_aic <- min(sel_table$AIC)
        # Zone of indifference (delta AIC < 2)
        top_tier <- sel_table[AIC <= (min_aic + 2)]
        
        # Parsimony Rule: Choose simplest model in the top tier
        setorder(top_tier, Complexity, AIC)
        
        best_idx <- top_tier$Index[1]
        best <- candidates[[best_idx]]
        
        if (nrow(top_tier) > 1) {
          ret$alternative_models_flag <- TRUE
        }
      }
    }
  }
  
  # Final check for model validity before emmeans
  if (inherits(best, "glmmTMB")) {
    bpue.r <- tryCatch(as.data.frame(emmeans(best, ~1, type="response", offset = log(1))),
                       error = function(e) NULL)
    if (!is.null(bpue.r)) {
      ret[, c("bpue", "lwr", "upr") := as.list(bpue.r[, c("response", "asymp.LCL", "asymp.UCL")])]
    }
  }
  
  # Format the model formula as a string
  if (inherits(best, "glmmTMB")) {
    ret[, model := paste(format(formula(best)), collapse = "")]
  } else {
    ret[, model := as.character(best)]
  }

	# In cases where no heterogeneity was detected, over-write best model by base model

  if (isFALSE(ret$base_model_heterogeneity)) {
    ret[, model := paste(format(formula(base_model)), collapse = "")]
  } else {}
									
    t_end <- Sys.time()
    t_elapsed <- difftime(t_end,t_start,units="mins")
    cat(sprintf("calc_bpue on %d rows completed in %.01f mins\n", nrow(needle), t_elapsed))

	bpue <- ret				   
    return(bpue)
    
}



