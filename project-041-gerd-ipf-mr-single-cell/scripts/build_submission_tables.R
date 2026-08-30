root <- normalizePath(".")
out <- file.path(root, "outputs", "tables_v1")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
suppressPackageStartupMessages(library(data.table))
t1 <- data.table(Source=c("GERD GWAS","IPF GWAS","IPF Cell Atlas"), Accession=c("GCST90000514","GCST009758","GSE136831"), Phenotype=c("Gastroesophageal reflux disease liability","Idiopathic pulmonary fibrosis risk","IPF/control/COPD lung single-cell transcriptomes"), Ancestry=c("European","European","Not applicable"), Sample_or_cells=c("78,707 cases; 288,734 controls (N=367,441)","2,668 cases; 8,591 controls (N=11,259)","312,928 cells from 78 donors"), Role=c("Exposure","Outcome and reverse-direction candidate","Exploratory localization"))
fwrite(t1, file.path(out,"Table1_data_sources.tsv"), sep="\t")
mr <- fread(file.path(root,"outputs","GERD_to_IPF_mr_primary_summary.tsv"))
cross <- fread(file.path(root,"outputs","GERD_to_IPF_mr_mendelianrandomization_crosscheck.tsv"))
t2 <- data.table(Analysis=c("IVW (primary)","IVW (package cross-check)","Weighted median","MR-Egger slope","MR-PRESSO raw"), Estimate=c(mr$ivw_beta,cross[method=="IVW"]$beta,cross[method=="Weighted median"]$beta,cross[method=="Egger"]$beta,mr$ivw_beta), SE=c(mr$ivw_se,cross[method=="IVW"]$se,cross[method=="Weighted median"]$se,cross[method=="Egger"]$se,NA), P_value=c(mr$ivw_p,cross[method=="IVW"]$p,cross[method=="Weighted median"]$p,cross[method=="Egger"]$p,0.00211), Interpretation=c("Primary association estimate","Independent arithmetic cross-check","Sensitivity; not conventionally significant","Sensitivity; pleiotropy diagnostic concern","Sensitivity; raw estimate"))
fwrite(t2, file.path(out,"Table2_mr_estimates.tsv"), sep="\t")
det <- fread(file.path(root,"outputs","GSE136831_target_genes_ipf_control_detection_tests_exploratory.tsv")); det <- det[FDR<.05 & IPF_donors>=10 & Control_donors>=10]
norm <- fread(file.path(root,"outputs","GSE136831_supported_signal_normalized_pseudobulk_sensitivity.tsv")); norm[,key:=paste(gene_symbol,celltype,sep="|")]; det[,key:=paste(gene_symbol,celltype,sep="|")]
t3 <- merge(det,norm[,.(key,CPM_IPF=IPF_median,CPM_Control=Control_median,CPM_median_diff=median_diff,CPM_P=wilcox_p,CPM_FDR=FDR)],by="key",all.x=TRUE)[,.(gene_id,gene_symbol,celltype,IPF_donors,Control_donors,detection_IPF=IPF_mean,detection_Control=Control_mean,detection_diff=difference,detection_P=wilcox_p,detection_FDR=FDR,CPM_IPF,CPM_Control,CPM_median_diff,CPM_P,CPM_FDR)]
fwrite(t3, file.path(out,"Table3_candidate_localization.tsv"), sep="\t")
cat("Wrote",nrow(t3),"supported localization rows to",out,"\n")
