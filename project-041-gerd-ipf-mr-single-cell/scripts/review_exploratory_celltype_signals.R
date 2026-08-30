library(data.table)

project <- normalizePath(".")
input <- file.path(project, "outputs", "GSE136831_target_genes_ipf_control_detection_tests_exploratory.tsv")
outdir <- file.path(project, "outputs")

t <- fread(input)
required <- c("gene_id", "gene_symbol", "celltype", "IPF_donors", "Control_donors",
              "IPF_mean", "Control_mean", "wilcox_p", "FDR")
stopifnot(all(required %in% names(t)))

t[, diff := IPF_mean - Control_mean]
t[, abs_diff := abs(diff)]
t[, support_gate := IPF_donors >= 10 & Control_donors >= 10]
t[, signal_gate := !is.na(FDR) & FDR < 0.05]
setorder(t, FDR, -abs_diff, gene_symbol, celltype)

fwrite(t, file.path(outdir, "GSE136831_exploratory_signal_review_all.tsv"), sep = "\t", na = "NA")

cand <- t[signal_gate & support_gate]
weak <- t[signal_gate & !support_gate]
fwrite(cand, file.path(outdir, "GSE136831_exploratory_signal_candidates_supported.tsv"), sep = "\t", na = "NA")
fwrite(weak, file.path(outdir, "GSE136831_exploratory_signal_candidates_weak_support.tsv"), sep = "\t", na = "NA")

audit <- c(
  "# Exploratory cell-type signal review v1",
  "",
  paste0("- Total rows: ", nrow(t)),
  paste0("- Rows with non-missing FDR: ", sum(!is.na(t$FDR))),
  paste0("- FDR<0.05 signals: ", sum(t$signal_gate)),
  paste0("- FDR<0.05 with >=10 donors in both IPF and Control: ", nrow(cand)),
  paste0("- FDR<0.05 with weaker donor support: ", nrow(weak)),
  paste0("- Cell types represented among supported candidates: ", uniqueN(cand$celltype)),
  "- Ranking includes absolute donor-level detection-prevalence difference.",
  "- Signals are exploratory and are not treated as confirmatory or causal evidence.",
  "- Supported candidates require sensitivity analysis, cell-type abundance review, and biological plausibility review before manuscript use.",
  "- Donor, not cell or library, remains the inferential unit."
)
writeLines(audit, file.path(outdir, "GSE136831_exploratory_signal_review_audit_v1.md"), useBytes = TRUE)

print(cand[, .(gene_symbol, celltype, IPF_donors, Control_donors, diff, wilcox_p, FDR)][1:min(20, .N)])
