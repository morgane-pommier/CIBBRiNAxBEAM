## ---------------------------
## Authors :  WGBYC ToR C members. Original scripts available at https://github.com/dlusseau/BEAM/tree/main
## Adapted for CIBBRiNA users by MP.
## ---------------------------
#
#'@description Assess zero bycatch probability when no incident/individual was observed, given the monitoring effort.
#'
#' @param bpue A data.table output from the BPUE estimation function (calc_bpue.r), containing at minimum columns: model, observed_effort, and the columns named in analysis_resolution. 
#' @param bycatch_rarity Numeric. Expected “true” probability of a bycatch incident for any given fishing operation. Default 0.001 (i.e. very rare, in 1000 fishing operations), as defined by WKPETSAMP (ICES, 2024) and agreed to in WKBEAM (ICES, 2026).   
#' @param analysis_resolution Character vector of column names defining the unit of analysis and the final aggregation level at which the assessment is produced. Must be the same as the resolution used to estimate BPUE (See calc_bpue.r).  
#' @return A data.table identical in structure to \code{bpue}, with bpue, lwr, and upr
#'   filled in for cases where zero bycatch is statistically supported (p <= 0.01).
#' @details A zero bpue is accepted for a given group when the probability of observing
#'   zero bycatch under \code{bycatch_rarity} given the observed effort is below 0.01,
#'   i.e. it is very unlikely that bycatch would have gone undetected. The upper confidence
#'   bound is set parametrically as 1.96 * sqrt(2 / (effort + 2)^2).
#'   Threshold values follow ICES WKBBEAM (Dec 2025), based on WKPETSAMP (ICES, 2024).

zeros_prob_vect <-function(n_ind=1L,effort,rarity=0.001,bootstrap=1000L) {
  mean(rbinom(bootstrap, size = effort * n_ind, prob = rarity) == 0L) #speed up the initial approach 
}

zero_assessments<-function(bpue,bycatch_rarity = 0.001, analysis_resolution ### Default thresholds decided at ICES WKBBEAM Dec 2025, based on WKPETSAMP report (ICES, 2024)
                           ) {
p <-0.01 #Threshold to accept / reject hypothesis. 0.01 chosen to be fairly certain that bycatch would have been detected and that therefore we can accept a bpue of 0.

zero.candidates <-bpue[bpue$model=="none",] #Only focus on cases where no bycatch data is available

zero.candidates[, effort := round(observed_effort, 0)]
zero.candidates <- zero.candidates[effort != 0]

zero.candidates$p_rarity<-unlist(lapply(zero.candidates$effort,function(x) zeros_prob_vect(n_ind=1L,effort=x,rarity=bycatch_rarity,bootstrap=1000L)))
zero.candidates$upper_parametric<-1.96 * sqrt(2 / (zero.candidates$effort+2)^2)

zc_ok <- unique(zero.candidates[p_rarity <= p, ..analysis_resolution])
zc_upper <- zero.candidates[, .(upper_parametric = upper_parametric), by = analysis_resolution]

bpue_with_zeros <- copy(bpue)

bpue_with_zeros[zc_ok, on = analysis_resolution, bpue := fifelse(is.na(bpue), 0, bpue)]
bpue_with_zeros[zc_ok, on = analysis_resolution, lwr := fifelse(is.na(bpue), 0, bpue)]
bpue_with_zeros[zc_upper, on = analysis_resolution, upr := i.upper_parametric]

return(bpue_with_zeros)
}
