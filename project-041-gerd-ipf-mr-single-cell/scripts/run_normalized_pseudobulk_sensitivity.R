library(data.table)

project <- normalizePath(".")
grid <- fread(file.path(project, "outputs", "GSE136831_target_genes_complete_donor_celltype_grid.tsv"))
meta <- fread(file.path(project, "inputs", "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"), quote='"', check.names=FALSE)
for (j in names(meta)) if (is.character(meta[[j]])) meta[[j]] <- gsub('"','',meta[[j]],fixed=TRUE)
libs <- fread(file.path(project, "outputs", "GSE136831_cell_library_sizes.tsv"))
stopifnot(nrow(meta) == nrow(libs), all(meta$CellBarcode_Identity == meta$CellBarcode_Identity))
meta[, cell_col := .I]
meta[, disease := Disease_Identity]
meta[, subject := Subject_Identity]
meta[, celltype := Subclass_Cell_Identity]

lib_group <- merge(meta[, .(cell_col, disease, subject, celltype)], libs, by="cell_col", all=FALSE)
lib_group <- lib_group[, .(library_umi = sum(library_umi), group_cells = .N), by=.(disease, subject, celltype)]
stopifnot(all(lib_group$library_umi > 0))

x <- merge(grid, lib_group, by=c("disease","subject","celltype"), all.x=TRUE, suffixes=c("_gene","_lib"))
stopifnot(nrow(x) == nrow(grid), all(!is.na(x$library_umi)))
x[, cpm := total_count / library_umi * 1e6]
x[, log_cpm1p := log1p(cpm)]
sig <- unique(fread(file.path(project, "outputs", "GSE136831_exploratory_signal_candidates_supported.tsv"),
                    select=c("gene_id","gene_symbol","celltype")))
x <- merge(x, sig, by=c("gene_id","gene_symbol","celltype"), all=FALSE)
x[, disease_group := fifelse(disease == "IPF", "IPF", fifelse(disease == "Control", "Control", NA_character_))]
x <- x[!is.na(disease_group)]

res <- x[, {
  a <- log_cpm1p[disease_group == "IPF"]; b <- log_cpm1p[disease_group == "Control"]
  p <- if (length(a) >= 10L && length(b) >= 10L && (length(unique(a)) > 1L || length(unique(b)) > 1L))
    wilcox.test(a,b,exact=FALSE)$p.value else NA_real_
  .(IPF_donors=length(a), Control_donors=length(b), IPF_median=median(a), Control_median=median(b), median_diff=median(a)-median(b), wilcox_p=p)
}, by=.(gene_id,gene_symbol,celltype)]
res[, FDR := p.adjust(wilcox_p, "BH")]
res[, abs_median_diff := abs(median_diff)]
setorder(res, FDR, -abs_median_diff)
fwrite(res, file.path(project, "outputs", "GSE136831_supported_signal_normalized_pseudobulk_sensitivity.tsv"), sep="\t", na="NA")
fwrite(lib_group, file.path(project, "outputs", "GSE136831_donor_celltype_library_umi.tsv"), sep="\t", na="NA")
writeLines(c(
  "# Normalized donor-level pseudobulk sensitivity audit v1", "",
  paste0("- Candidate signals evaluated: ", nrow(res)),
  paste0("- Non-missing normalized P values: ", sum(!is.na(res$wilcox_p))),
  paste0("- Normalized FDR<0.05: ", sum(!is.na(res$FDR) & res$FDR < 0.05)),
  "- Normalization: candidate-gene pseudobulk counts divided by total cell-library UMI summed within each donor-cell type (CPM), followed by log1p.",
  "- The full matrix was streamed to obtain cell-library UMI; no cell-level inferential test was performed.",
  "- This is a donor-level sensitivity analysis; it is not a DESeq2/edgeR model using genome-wide counts and should not be described as definitive differential expression.",
  "- Concordance with detection prevalence is supportive only and remains subject to cell-composition and donor-level confounding."
), file.path(project, "outputs", "GSE136831_normalized_pseudobulk_sensitivity_audit_v1.md"), useBytes=TRUE)
print(res)
