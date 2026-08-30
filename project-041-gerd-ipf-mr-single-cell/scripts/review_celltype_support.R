library(data.table)

project <- normalizePath(".")
grid <- fread(file.path(project, "outputs", "GSE136831_target_genes_complete_donor_celltype_grid.tsv"))
cand <- fread(file.path(project, "outputs", "GSE136831_exploratory_signal_candidates_supported.tsv"))

subject_level <- unique(grid[, .(disease, subject, celltype, group_cells = group_cells.x)])
subject_level[, group_cells := as.numeric(group_cells)]
subject_level[, disease_group := fifelse(disease == "IPF", "IPF", fifelse(disease == "Control", "Control", disease))]
coverage <- subject_level[disease_group %in% c("IPF", "Control"),
  .(donors = uniqueN(subject), median_cells = median(group_cells),
    q1_cells = quantile(group_cells, 0.25, names = FALSE),
    q3_cells = quantile(group_cells, 0.75, names = FALSE),
    min_cells = min(group_cells), max_cells = max(group_cells)),
  by = .(celltype, disease_group)]
setorder(coverage, celltype, disease_group)
fwrite(coverage, file.path(project, "outputs", "GSE136831_celltype_donor_coverage_summary.tsv"), sep = "\t", na = "NA")

candidate_support <- merge(cand, coverage, by = "celltype", all.x = TRUE, allow.cartesian = TRUE)
candidate_support[, coverage_flag := (disease_group == "IPF" & donors == IPF_donors) |
                                  (disease_group == "Control" & donors == Control_donors)]
fwrite(candidate_support, file.path(project, "outputs", "GSE136831_supported_signal_celltype_context.tsv"), sep = "\t", na = "NA")

writeLines(c(
  "# Cell-type donor coverage audit v1", "",
  paste0("- Candidate signals reviewed: ", nrow(cand)),
  paste0("- Unique candidate cell types: ", uniqueN(cand$celltype)),
  paste0("- Coverage rows (IPF/Control): ", nrow(coverage)),
  "- Coverage and cell-count summaries are descriptive and do not adjust exploratory P values.",
  "- Candidate signals require composition and donor-level count sensitivity checks before manuscript use.",
  "- Donor, not cell or library, remains the inferential unit."
), file.path(project, "outputs", "GSE136831_celltype_donor_coverage_audit_v1.md"), useBytes = TRUE)

print(coverage[celltype %in% unique(cand$celltype)])
