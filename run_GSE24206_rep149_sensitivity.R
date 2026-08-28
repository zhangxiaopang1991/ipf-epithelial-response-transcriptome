sc <- read.csv("outputs/mechanism_module_scores.csv", check.names=FALSE, stringsAsFactors=FALSE)
module <- "fibrosis_epithelial_response"
d <- sc[sc$dataset=="GSE24206" & sc$group %in% c("early_IPF","control"),]
drop_ids <- c("none", "GSM595426", "GSM595427")
rows <- lapply(drop_ids, function(id) {
  x <- if (id=="none") d else d[d$sample_id != id,]
  a <- x[[module]][x$group=="early_IPF"]
  b <- x[[module]][x$group=="control"]
  w <- wilcox.test(a,b,exact=FALSE)
  data.frame(excluded_sample=ifelse(id=="none",NA_character_,id), n_early=length(a), n_control=length(b),
    median_early=median(a), median_control=median(b), median_delta=median(a)-median(b), p=w$p.value,
    direction=ifelse(median(a)>median(b),"early_higher","early_lower"), stringsAsFactors=FALSE)
})
r <- do.call(rbind, rows)
stopifnot(all(r$direction=="early_higher"))
write.csv(r,"outputs/GSE24206_rep149_sample_level_sensitivity.csv",row.names=FALSE,fileEncoding="UTF-8")
capture.output({
  cat("GSE24206 sample-level sensitivity for duplicate donor label rep149\n")
  cat("GSM595426 = upper-lobe biopsy; GSM595427 = lower-lobe biopsy; both are early IPF, male, age 58.\n")
  cat("This sensitivity treats one sample at a time as removable. It does not recover a patient-level analysis.\n\n")
  print(r)
},file="outputs/GSE24206_rep149_sample_level_sensitivity_report.txt")
