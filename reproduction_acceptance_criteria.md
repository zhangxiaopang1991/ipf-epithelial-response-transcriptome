# Reproduction acceptance criteria

Run the scripts from a clean working directory containing only `scripts/`, empty `outputs/`, and the documented public GEO inputs.

The reproduction passes when all conditions below hold:

1. `verify_discovery_validation_package.R` prints `DISCOVERY_VALIDATION_ASSERTIONS_PASSED`.
2. `main_results_package_audit.txt` begins with `Main results package audit passed.`
3. The discovery-direction comparisons have sample counts 119/50, 23/15, and 8/6, respectively, with positive module-score median differences.
4. The GSE32537 adjusted, FVC, and DLCO models have n=169, n=117, and n=99, respectively.
5. The GSE32537 primary-module coefficients are approximately 1.103, -13.314, and -15.272, respectively, at the precision reported in the manuscript.
6. The primary-module robust coefficient is approximately 0.959. Differences below 1e-9 in floating-point estimates or P values are acceptable because robust regression can vary at machine precision across otherwise identical executions.
7. GO includes `extracellular matrix organization` with 17/212 genes and FDR approximately 0.0146. Gene-list order within the `geneID` field is not an inferential result and may differ.
8. `run_limma_gsea_GSE32537.R` always writes `limma_GSE32537_IPF_UIP_vs_control_adjusted.csv`. Hallmark GSEA is optional; a missing external gene-set download must be reported as skipped and must not block GO.

No criterion permits an inference about reflux, laryngopharyngeal reflux, microaspiration, laryngoscopy, causality, prognosis, prediction, or longitudinal progression.
