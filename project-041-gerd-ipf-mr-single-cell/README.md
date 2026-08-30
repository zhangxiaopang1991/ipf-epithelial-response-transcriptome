# Project 041 code release candidate

This directory contains the reproducible code for the GERD-IPF MR plus IPF lung single-cell localization study.

## Scope

- Two-sample MR using GCST90000514 and GCST009758.
- Exploratory donor-level localization in GSE136831.
- Library-size-normalized candidate-gene CPM sensitivity analysis.
- Reproducible figure-generation and table-generation scripts.

## Important boundaries

The single-cell module is exploratory localization, not causal mediation or confirmatory differential expression. The reverse MR direction is limited by instrument availability. No participant-level data are included.

## Run order

Run commands from this directory. Public input files and the project-local `.r_libs` directory are not included in the repository.

1. `run_ld_clump.R`
2. `annotate_mr_snps_ensembl.ps1` and `build_candidate_gene_list.ps1`
3. `build_pre_mr_input_table.ps1`
4. `run_gerd_to_ipf_mr_baseR.R`, `crosscheck_mendelianrandomization.R` and `run_mrpresso_gerd_to_ipf.R`
5. `audit_gse136831_metadata.R`, `check_candidate_gene_index.R` and `extract_matched_genes_stream.js`
6. `aggregate_target_genes_by_donor_celltype.R`, `complete_and_summarize_target_genes.R`
7. `review_exploratory_celltype_signals.R`, `review_celltype_support.R`
8. `run_expression_intensity_sensitivity.R`, `compute_cell_library_sizes.js`, `run_normalized_pseudobulk_sensitivity.R`
9. `make_submission_figures.R` and `build_submission_tables.R`

## Before release maintenance

1. Run the scripts from a clean environment and record the exact session information.
2. Check that all scripts use relative paths and contain no local credentials.
3. Compare regenerated tables and figures with the audited outputs.
4. Keep the exploratory and reverse-MR limitations in any derived report.
