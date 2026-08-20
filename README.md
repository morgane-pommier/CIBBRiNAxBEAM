# BEAM-CIBBRiNA

This repository compiles functions from the CIBBRiNAxBEAM bycatch estimation toolbox. The tools are implemented in the R programming language. The core functions of ICES WGBYC BEAM, used to estimate bycatch per unit of effort (BPUE) and total bycatch, are adapted to facilitate their application to a wider and more diverse set of bycatch data, such as those collected in CIBBRiNA case studies: 

- **calc_bpue.r**: Assesses whether heterogeneity in BPUE (Bycatch per Unit of Effort) occurs within the fleet, selects the best modelling option, and estimates population-level BPUEs.
- **zero_assessment.r**: For cases where no bycatch was reported, assesses whether a BPUE of 0 can be statistically supported, given the monitoring effort available.
- **calc_total.r**: Uses the best model identified in calc_bpue.r to raise BPUE to total fishing effort and compute total bycatch, where possible.
- **calc_partia.r**: In cases where total bycatch cannot be estimated due to variability within the fleet, but incomplete monitoring coverage, calculates a partial bycatch estimate for the part of the fishery in which BPUE is known.
- **reliability_estimation.r**: Assesses the reliability of total bycatch estimates by combining a measure of BPUE robustness (i.e. how stable the estimate is when individual data points are iteratively from model fitting) , with the breadth of the confidence interval around the total bycatch estimate.

For a tutorial on how to use the toolbox, please visit: https://morgane-pommier.github.io/CIBBRiNAxBEAM/. This implementation example uses a simulated dataset available in the /data folder, and all the code is provided in the Tutorial.Rmd file.

The report detailing those functions is available at **[INSERT REF TO DELIVERABLE WHEN AVAILABLE]**

## BEAM : Bycatch Estimation and Assessment Matrix, by ICES WGBYC.

BEAM is the estimation and assessment method used by ICES Working Group on Bycatch of Endangered, Threatened and Protected species (WGBYC, [ICES, 2026](https://ices-library.figshare.com/articles/report/Benchmark_Workshop_on_Bycatch_Evaluation_and_Assessment_Matrix_WKBBEAM_/32922134?file=66384137))

The original idea of BEAM was developed by Marjorie Lyssikatos and David Lusseau. The code was developed by David Lusseau, André Moan, Paula Gutiérrez Muñoz, Morgane Pommier, Henrik Pärn, Torbjörn Säterberg, Kim Magnus Bærum and Amaia Astarloa Diaz. WGBYC ToR C members contributed significantly to non-coding development and beta testing, particularly Marjorie Lyssikatos, Ailbhe Kavanagh, Caterina Fortuna, Guðjón Már Sigurðsson, Katja Ringdahl, Sara Königson, Ruth Fernandez, Allen Kingston, Lotte Kindt-Larsen and Carlos Pinto. The efforts to get BEAM to ICES Transparent Assessment Framework (TAF) are led by André Moan.
The original BEAM code, used for ICES annual assessment by WGBYC and adapted here for CIBBRiNA, is available at: https://github.com/dlusseau/BEAM. 

Further developments under CIBBRiNA WP6 were led by Morgane Pommier, Ailbhe Kavanagh and David Lusseau.

## Advice on data preparation and quality control

To use the CIBBRiNAxBEAM toolbox, users should prepare two datasets: one combining monitoring effort and bycatch records, to estimate BPUE (Bycatch per Unit Effort), and one describing total fishing effort to scale BPUE up to total bycatch in the fishery. Both datasets should be aggregated at the same resolution. 

<img width="5795" height="2575" alt="Data structure" src="https://github.com/user-attachments/assets/8ead43da-6743-4d9a-8e9e-c675ff2be3cf" />

Aggregation levels (Figure 1) should capture the scale at which the assessment needs to be produced (e.g. species x gear x region), while also accounting for finer grouping levels that may influence bycatch rate (e.g. year x subregion x vessel ID or vessel size x season x target fishery etc). Additional variables describing fishing operations, such as soak time or net length, may also be included, alongside the main unit of effort (e.g. number of days, trips, hooks etc) and corresponding observed bycatch. 
Users may wish to apply quality standards and data filtering steps before implementing the BEAM process. For instance, in an ICES assessment, logbook data and port-observer data are excluded, as are records made while the protocol was focused on a different taxon than the one of the species of interest.

Another recommended quality control step is to verify that in each substratum (i.e. the finest level of aggregation available), the monitoring effort is lower than (or equal to) the total fishing effort. If it’s higher, this may indicate that not all fishing effort was reported or that monitoring data needs to be corrected. While the procedure can still be performed and BPUE values won’t be affected, it is worth acknowledging potential issues when discussing the results, as it may affect the total bycatch estimate (underestimated total fishing effort will likely result in underestimated bycatch).


