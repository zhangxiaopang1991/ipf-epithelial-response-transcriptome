library(data.table)
project <- normalizePath(".")
agg <- fread(file.path(project,"outputs","GSE136831_target_genes_donor_celltype_aggregation.tsv"))
groups <- fread(file.path(project,"outputs","GSE136831_donor_celltype_group_sizes.tsv"))
genes <- unique(agg[, .(gene_id, gene_symbol)])
grid <- as.data.table(merge(as.data.frame(genes), as.data.frame(groups), by=NULL))
grid <- merge(grid, agg, by=c("gene_id","gene_symbol","disease","subject","celltype"), all.x=TRUE)
grid[is.na(total_count), `:=`(total_count=0, detected_cells=0, detection_fraction=0, mean_count_per_cell=0)]
setorder(grid, gene_symbol, disease, subject, celltype)
fwrite(grid, file.path(project,"outputs","GSE136831_target_genes_complete_donor_celltype_grid.tsv"), sep="\t")
summary <- grid[, .(donors=.N, donors_detected=sum(detection_fraction>0), mean_detection_fraction=mean(detection_fraction), median_detection_fraction=median(detection_fraction), mean_count_per_cell=mean(mean_count_per_cell)), by=.(gene_id,gene_symbol,disease,celltype)]
summary[, detection_prevalence := donors_detected / donors]
fwrite(summary, file.path(project,"outputs","GSE136831_target_genes_detection_summary.tsv"), sep="\t")
donor_prev <- grid[, .(detection_prevalence=mean(detection_fraction)), by=.(gene_id,gene_symbol,disease,subject,celltype)]
ipf_ctrl <- donor_prev[disease %in% c("IPF","Control")]
tests <- ipf_ctrl[, { x <- detection_prevalence[disease=="IPF"]; y <- detection_prevalence[disease=="Control"]; n1 <- length(x); n0 <- length(y); p <- if(n1>=5 && n0>=5) suppressWarnings(wilcox.test(x,y,exact=FALSE)$p.value) else NA_real_; .(IPF_donors=n1,Control_donors=n0,IPF_mean=if(n1)mean(x) else NA_real_,Control_mean=if(n0)mean(y) else NA_real_,wilcox_p=p) }, by=.(gene_id,gene_symbol,celltype)]
tests[, FDR := p.adjust(wilcox_p, method="BH")]
fwrite(tests, file.path(project,"outputs","GSE136831_target_genes_ipf_control_detection_tests_exploratory.tsv"), sep="\t")
writeLines(c("# GSE136831 complete-grid and detection summary audit v1","",paste0("- Complete grid rows: ",nrow(grid)),paste0("- Genes: ",uniqueN(grid$gene_symbol)),paste0("- Donor-celltype groups: ",nrow(groups)),paste0("- Missing aggregate values after zero fill: ",sum(is.na(grid$total_count)|is.na(grid$detection_fraction))),"- Exploratory IPF-vs-Control tests use donor-level detection prevalence and require at least 5 donors per group; BH FDR is descriptive.","- No cell-level inferential test was performed."), file.path(project,"outputs","GSE136831_complete_grid_audit_v1.md"))
print(grid[, .(rows=.N, zero_filled=sum(total_count==0)), by=disease])
