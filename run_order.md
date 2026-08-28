# Reproduction run order

Run from the repository root, after placing the public GEO matrix and GPL annotation files in the documented `inputs/` paths. Do not place raw files in a public code repository unless the repository and GEO terms permit it.

1. `Rscript scripts/audit_geo_metadata.R`
2. `Rscript scripts/run_mechanism_pilot.R`
3. `Rscript scripts/run_main_module_robustness.R`
4. `Rscript scripts/run_adjusted_models.R`
5. `Rscript scripts/run_robust_regression_sensitivity.R`
6. `Rscript scripts/run_gsva_cell_sensitivity.R`
7. `Rscript scripts/run_limma_gsea_GSE32537.R`; adjusted limma output is required by GO, while Hallmark GSEA is optional and skips cleanly if its gene-set dependency is unavailable
8. `Rscript scripts/run_go_ora_GSE32537.R`
9. `Rscript scripts/run_GSE24206_rep149_sensitivity.R`
10. `Rscript scripts/run_stratified_case_definitions.R`
11. `Rscript scripts/run_validation_summary.R`
12. `Rscript scripts/build_discovery_validation_package.R`
13. `Rscript scripts/verify_discovery_validation_package.R`
14. `Rscript scripts/build_main_results_package.R`
15. `Rscript scripts/make_pilot_figure.R` (legacy combined pilot figure)
16. `Rscript scripts/make_submission_figures.R` (Figure 1-3 and Supplementary Figure S1)

The manuscript's primary claims require successful completion of steps 2, 3, 4, 5, 6, 8, 9, 12, 13, and 14. Step 7 is optional because Hallmark GSEA was not reported. `run_cell_proxy_sensitivity.R` is retained as an exploratory legacy sensitivity script but is not part of the reported primary reproduction order.
