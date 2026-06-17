# BEAM-CIBBRiNA

This repository compiles functions from the CIBBRiNAxBEAM bycatch estimation toolbox. The core functions of ICES WGBYC BEAM (https://github.com/dlusseau/BEAM), used to estimate bycatch per unit of effort (BPUE) and total bycatch, were further developed to facilitate their application to a wider and more diverse set of bycatch data, such as those collected in CIBBRiNA case studies. 

- calc_bpue.r: Assesses whether heterogeneity in BPUE occurs within the fleet, selects the best modelling option in consequence, and estimates population-level BPUEs.
- zero_assessment.r: For cases where no bycatch was reported, assesses whether a BPUE of 0 can be statistically supported, given the monitoring effort available.
- calc_total.r: Uses the best model identified in calc_bpue.r to raise BPUE to total fishing effort and compute total bycatch, where possible.
- calc_partia.r: In cases where total bycatch could not be estimated due to variability in BPUE across fleet segments but incomplete sampling across segments, calculates a partial bycatch estimate for the part of the fishery in which BPUE is known.
- reliability_estimation.r: Assesses the reliability of total bycatch estimates by combining a measure of BPUE robustness to exclusion of individual observations (based on RMSE), with the breadth of the confidence interval around the total bycatch estimate.

The report detailing those functions is available at [INSERT REF TO DELIVERABLE]

# Advice on data preparation and quality control

To use these functions, users should prepare two datasets: one combining monitoring effort and bycatch records to estimate BPUE, and one describing total fishing effort to scale BPUE up to total bycatch. Both datasets should be aggregated at the same resolution. 

<img width="982" height="423" alt="image" src="https://github.com/user-attachments/assets/fa89a0a2-d0e4-4453-8869-f32230a580f7" />

Aggregation levels (Figure 1) should capture the scale at which the assessment is to be produced (e.g. species x gear x region), while also accounting for finer grouping levels that may influence bycatch rate (e.g. year x subregion x vessel ID or vessel size x season x target fishery etc). Additional variables describing fishing operations, such as soak time or net length, may also be included, alongside the main unit of effort (e.g. number of days, trips, hooks etc) and corresponding observed bycatch. 
Users may wish to apply quality standards and data filtering steps before implementing the BEAM process. For instance, in an ICES assessment, logbook data and port-observer data are excluded, as are records made while the protocol was focused on a different taxon than the one of the species of interest.

Another recommended quality control step is to verify that in each substratum (i.e. the finest level of aggregation available), the monitoring effort is lower than (or equal to) the total fishing effort. If it’s higher, this may indicate that not all fishing effort was reported or that monitoring data needs to be corrected. While the procedure can still be performed and BPUE values won’t be affected, it is worth acknowledging potential issues when discussing the results, as it may affect the total bycatch estimate (underestimated total fishing effort will likely result in underestimated bycatch).


