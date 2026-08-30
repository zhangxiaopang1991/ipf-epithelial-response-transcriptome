library(data.table)

project <- "workspace-code-sop-project-agent-md/projects/project-041-ipf-gerd-mr-single-cell"
grid <- fread(file.path(project, "outputs", "GSE136831_target_genes_complete_donor_celltype_grid.tsv"))
sig <- fread(file.path(project, "outputs", "GSE136831_exploratory_signal_candidates_supported.tsv"),
             select = c("gene_id", "gene_symbol", "celltype"))
sig <- unique(sig)

x <- merge(grid, sig, by = c("gene_id", "gene_symbol", "celltype"), all = FALSE)
x[, disease_group := fifelse(disease == "IPF", "IPF", fifelse(disease == "Control", "Control", NA_character_))]
x <- x[!is.na(disease_group)]
x[, intensity_log1p := log1p(as.numeric(mean_count_per_cell))]

res <- x[, {
  ipf <- intensity_log1p[disease_group == "IPF"]
  ctrl <- intensity_log1p[disease_group == "Control"]
  p <- if (length(ipf) >= 10L && length(ctrl) >= 10L && (length(unique(ipf)) > 1L || length(unique(ctrl)) > 1L))
    wilcox.test(ipf, ctrl, exact = FALSE)$p.value else NA_real_
  .(IPF_donors = length(ipf), Control_donors = length(ctrl),
    IPF_median = median(ipf), Control_median = median(ctrl),
    median_diff = median(ipf) - median(ctrl), wilcox_p = p)
}, by = .(gene_id, gene_symbol, celltype)]
res[, FDR := p.adjust(wilcox_p, method = "BH")]
res[, abs_median_diff := abs(median_diff)]
setorder(res, FDR, -abs_median_diff)

fwrite(res, file.path(project, "outputs", "GSE136831_supported_signal_expression_intensity_sensitivity.tsv"), sep = "\t", na = "NA")
writeLines(c(
  "# Expression-intensity sensitivity audit v1", "",
  paste0("- Supported detection-prevalence signals evaluated: ", nrow(res)),
  paste0("- Non-missing intensity P values: ", sum(!is.na(res$wilcox_p))),
  paste0("- Intensity FDR<0.05: ", sum(!is.na(res$FDR) & res$FDR < 0.05)),
  "- Metric: donor-level log1p(mean raw count per cell) within each cell type.",
  "- This is an unnormalized descriptive sensitivity analysis and is not a substitute for donor-aware normalized pseudobulk modeling.",
  "- Concordance with detection prevalence is supportive only; discordance should be reported as exploratory uncertainty."
), file.path(project, "outputs", "GSE136831_expression_intensity_sensitivity_audit_v1.md"), useBytes = TRUE)
print(res)
