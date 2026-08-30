library(data.table)
project <- normalizePath(".")
meta <- fread(file.path(project,"inputs","GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"), quote='"', check.names=FALSE)
for (j in names(meta)) if (is.character(meta[[j]])) meta[[j]] <- gsub('"','',meta[[j]],fixed=TRUE)
bc_con <- gzfile(file.path(project,"inputs","GSE136831_AllCells.cellBarcodes.txt.gz"),"rt")
barcodes <- readLines(bc_con,warn=FALSE)
close(bc_con)
barcodes <- gsub('"','',barcodes,fixed=TRUE)
stopifnot(nrow(meta) == length(barcodes), nrow(meta) == 312928)
stopifnot(identical(as.character(meta$CellBarcode_Identity), barcodes))
meta[, cell_col := .I]
meta[, disease := Disease_Identity]
meta[, subject := Subject_Identity]
meta[, celltype := Subclass_Cell_Identity]
group_sizes <- meta[, .(group_cells=.N), by=.(disease, subject, celltype)]
sparse <- fread(file.path(project,"outputs","GSE136831_matched_genes_sparse.tsv"))
merged <- merge(sparse, meta[,.(cell_col, disease, subject, celltype)], by="cell_col", all.x=TRUE)
stopifnot(nrow(merged) == nrow(sparse), all(!is.na(merged$subject)))
agg <- merged[, .(total_count=sum(count), detected_cells=sum(count > 0)), by=.(gene_id, gene_symbol, disease, subject, celltype)]
agg <- merge(agg, group_sizes, by=c("disease","subject","celltype"), all.x=TRUE)
agg[, detection_fraction := detected_cells / group_cells]
agg[, mean_count_per_cell := total_count / group_cells]
setorder(agg, gene_symbol, disease, subject, celltype)
fwrite(agg, file.path(project,"outputs","GSE136831_target_genes_donor_celltype_aggregation.tsv"), sep="\t")
fwrite(group_sizes, file.path(project,"outputs","GSE136831_donor_celltype_group_sizes.tsv"), sep="\t")
summary <- agg[, .(genes=uniqueN(gene_symbol), donor_celltype_rows=.N, donors=uniqueN(subject), celltypes=uniqueN(celltype)), by=disease]
fwrite(summary, file.path(project,"outputs","GSE136831_target_genes_aggregation_summary.tsv"), sep="\t")
writeLines(c("# GSE136831 donor-celltype aggregation audit v1","",paste0("- Sparse rows aggregated: ",nrow(agg)),paste0("- Genes represented: ",uniqueN(agg$gene_symbol),"/54"),paste0("- Donor-celltype groups represented: ",nrow(group_sizes)),paste0("- Missing denominator group sizes after merge: ",sum(is.na(agg$group_cells))),"- Cell columns were matched to metadata barcodes by exact order and identity assertion."), file.path(project,"outputs","GSE136831_target_genes_aggregation_audit_v1.md"))
print(summary)
