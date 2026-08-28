out <- normalizePath("outputs", mustWork=TRUE)
adj <- read.csv(file.path(out,"adjusted_linear_models_GSE32537.csv"))
rob <- read.csv(file.path(out,"robust_lmrob_no_tissue_source_GSE32537.csv"))
strat <- read.csv(file.path(out,"stratified_case_definition_tests.csv"))
loo <- read.csv(file.path(out,"main_module_leave_one_out_and_gsva.csv"))
main <- subset(adj, analysis=="IPF_UIP_vs_control_adjusted" & outcome=="fibrosis_epithelial_response")
fvc <- subset(adj, outcome=="FVC_percent" & grepl("fibrosis_epithelial_response$",analysis))
dlco <- subset(adj, outcome=="DLCO_percent" & grepl("fibrosis_epithelial_response$",analysis))
robm <- subset(rob, outcome=="fibrosis_epithelial_response")
stratm <- subset(strat, module=="fibrosis_epithelial_response")
full <- subset(loo, variant=="full")
package <- rbind(
 data.frame(domain="GSE32537_adjusted_case_control", comparison="IPF_UIP_vs_control", n=main$n, estimate=main$beta, p=main$p, fdr=main$fdr_within_all_models, evidence="adjusted linear model"),
 data.frame(domain="GSE32537_robust_sensitivity", comparison="IPF_UIP_vs_control", n=robm$n, estimate=robm$beta, p=robm$p, fdr=robm$fdr, evidence="lmrob without tissue source"),
 data.frame(domain="GSE32537_lung_function", comparison="FVC_perc", n=fvc$n, estimate=fvc$beta, p=fvc$p, fdr=fvc$fdr_within_all_models, evidence="cross-sectional adjusted linear model"),
 data.frame(domain="GSE32537_lung_function", comparison="DLCO_perc", n=dlco$n, estimate=dlco$beta, p=dlco$p, fdr=dlco$fdr_within_all_models, evidence="cross-sectional adjusted linear model"),
 data.frame(domain="cross_dataset_score_robustness", comparison=paste(full$dataset,full$method,sep="_"), n=full$n_case+full$n_control, estimate=full$delta, p=full$p, fdr=NA, evidence="predefined module score"),
 data.frame(domain="stratified_case_definition", comparison=paste(stratm$dataset,stratm$case_group,"vs",stratm$reference_group), n=stratm$n_case+stratm$n_reference, estimate=stratm$delta, p=stratm$p, fdr=stratm$fdr_within_contrast, evidence="exploratory Wilcoxon")
)
write.csv(package,file.path(out,"main_results_package.csv"),row.names=FALSE,fileEncoding="UTF-8")
stopifnot(main$n==169, robm$n==169, fvc$n==117, dlco$n==99, all(full$delta>0), all(stratm$delta[stratm$case_group %in% c("UIP","early_IPF","advanced_IPF","IPF_UIP")]>0))
capture.output({cat("Main results package audit passed.\n");print(package)},file=file.path(out,"main_results_package_audit.txt"))
