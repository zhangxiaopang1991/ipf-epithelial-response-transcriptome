set.seed(20260827)
sc <- read.csv("outputs/mechanism_module_scores.csv", check.names=FALSE, stringsAsFactors=FALSE)
module <- "fibrosis_epithelial_response"

boot_delta <- function(a, b, B=5000) {
  f <- function() median(sample(a, length(a), replace=TRUE)) - median(sample(b, length(b), replace=TRUE))
  x <- replicate(B, f())
  as.numeric(quantile(x, c(.025,.975), names=FALSE, type=6))
}
rank_biserial <- function(a,b) {
  z <- outer(a,b,"-")
  (sum(z>0) + 0.5*sum(z==0) - 0.5*sum(z<0)) / (length(a)*length(b))
}
contrast <- data.frame(
  dataset=c("GSE32537","GSE10667","GSE24206"),
  role=c("discovery","directional_validation","directional_validation"),
  case_group=c("IPF_UIP","UIP","early_IPF"),
  reference_group=c("control","control","control"),
  stringsAsFactors=FALSE)
rows <- vector("list", nrow(contrast))
for (i in seq_len(nrow(contrast))) {
  q <- contrast[i,]
  d <- sc[sc$dataset==q$dataset & sc$group %in% c(q$case_group,q$reference_group),]
  a <- d[[module]][d$group==q$case_group]; b <- d[[module]][d$group==q$reference_group]
  ci <- boot_delta(a,b)
  wt <- wilcox.test(a,b, exact=FALSE)
  rows[[i]] <- data.frame(dataset=q$dataset, role=q$role, case_group=q$case_group,
    reference_group=q$reference_group, n_case=length(a), n_reference=length(b),
    median_case=median(a), median_reference=median(b), median_delta=median(a)-median(b),
    bootstrap_ci_low=ci[1], bootstrap_ci_high=ci[2], rank_biserial=rank_biserial(a,b),
    wilcoxon_p=wt$p.value, stringsAsFactors=FALSE)
}
r <- do.call(rbind, rows)
write.csv(r, "outputs/discovery_validation_main_module.csv", row.names=FALSE, fileEncoding="UTF-8")
capture.output({
  cat("Discovery-validation package for prespecified fibrosis/epithelial-response module\n")
  cat("Score: platform-wise gene-wise z-score mean; bootstrap B=5000; seed=20260827\n\n")
  print(r)
  cat("\nInterpretation boundary:\n")
  cat("GSE32537 is the discovery/clinical-association cohort. GSE10667 and GSE24206 test only directional replication of the fixed module in prespecified case-control contrasts; they are not independent clinical prediction validation cohorts.\n")
}, file="outputs/discovery_validation_main_module_report.txt")
