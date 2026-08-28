library(GEOquery)
library(limma)

base <- normalizePath(file.path("inputs", "geo_series_matrix"), mustWork = TRUE)
out <- normalizePath("outputs", mustWork = TRUE)
ann_dir <- normalizePath("inputs", mustWork = TRUE)

modules <- list(
  epithelial_barrier = c("EPCAM","CDH1","OCLN","CLDN18","CLDN4","KRT8","KRT18"),
  mucociliary_clearance = c("FOXJ1","PIFO","DNAH5","DNAI1","MUC1","MUC5B","SCGB1A1"),
  injury_inflammation = c("IL1B","IL6","CXCL8","NFKB1","TNF","TLR4"),
  oxidative_stress = c("NFE2L2","HMOX1","SOD2","GPX2","NOX4"),
  fibrosis_epithelial_response = c("TGFB1","COL1A1","COL3A1","SERPINE1","VIM","MMP7","KRT19")
)

read_annotation <- function(gpl) {
  f <- file.path(ann_dir, paste0(gpl, ".annot.gz"))
  con <- gzfile(f, "rt")
  repeat { z <- readLines(con, n = 1, warn = FALSE); if (!length(z)) stop("annotation table not found: ", gpl); if (z == "!platform_table_begin") break }
  lines <- readLines(con, warn = FALSE)
  close(con)
  if (!length(lines)) stop("empty annotation table: ", gpl)
  hdr <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  x <- read.delim(text = paste(lines[-1], collapse = "\n"), header = FALSE, sep = "\t", quote = '"', fill = TRUE, check.names = FALSE, col.names = hdr, comment.char = "")
  x[, c("ID", "Gene symbol")]
}

read_eset <- function(gse) getGEO(filename = file.path(base, paste0(gse, "_series_matrix.txt.gz")), GSEMatrix = TRUE, getGPL = FALSE)

parse_group <- function(pd, gse) {
  char_cols <- grep("^characteristics", names(pd), value = TRUE)
  tx <- tolower(apply(pd[, c("title", "source_name_ch1", char_cols), drop = FALSE], 1, paste, collapse = " | "))
  if (gse == "GSE24206") return(factor(ifelse(grepl("advanced|late|transplant", tx), "advanced_IPF", ifelse(grepl("early", tx), "early_IPF", ifelse(grepl("healthy|normal|control", tx), "control", "other")))))
  if (gse == "GSE10667") return(factor(ifelse(grepl("aex|acute", tx), "AEx", ifelse(grepl("uip|ipf", tx), "UIP", ifelse(grepl("control|normal", tx), "control", "other")))))
  factor(ifelse(grepl("ipf|uip", tx), "IPF_UIP", ifelse(grepl("control|normal", tx), "control", "other")))
}

score_dataset <- function(gse) {
  es <- read_eset(gse); e <- exprs(es); pd <- pData(es); gpl <- as.character(pd$platform_id[1])
  ann <- read_annotation(gpl); ann <- ann[match(rownames(e), ann$ID), ]
  sym <- toupper(trimws(sub("///.*$", "", as.character(ann$`Gene symbol`))))
  keep <- !is.na(sym) & sym != ""
  e <- e[keep, , drop = FALSE]; sym <- sym[keep]
  e <- avereps(e, ID = sym)
  present <- unique(sym)
  scores <- data.frame(dataset = gse, sample_id = colnames(e), group = parse_group(pd, gse), stringsAsFactors = FALSE)
  for (nm in names(modules)) {
    genes <- intersect(modules[[nm]], rownames(e))
    if (!length(genes)) { scores[[nm]] <- NA_real_; next }
    z <- t(scale(t(e[genes, , drop = FALSE])))
    scores[[nm]] <- colMeans(z, na.rm = TRUE)
  }
  if (gse == "GSE32537") {
    cc <- tolower(apply(pd[, grep("^characteristics", names(pd)), drop = FALSE], 1, paste, collapse = " | "))
    get_num <- function(pattern) {
      m <- regexec(paste0(pattern, "[^0-9]*([0-9]+(\\.[0-9]+)?)"), cc)
      vapply(regmatches(cc, m), function(z) if (length(z) >= 2) as.numeric(z[2]) else NA_real_, numeric(1))
    }
    scores$FVC_percent <- get_num("fvc")
    scores$DLCO_percent <- get_num("dlco")
  }
  list(scores = scores, genes = data.frame(dataset = gse, platform = gpl, module = names(modules), n_genes = vapply(modules, function(x) length(intersect(x, present)), integer(1))))
}

res <- lapply(c("GSE10667", "GSE24206", "GSE32537"), score_dataset)
score_list <- lapply(res, `[[`, "scores")
all_cols <- unique(unlist(lapply(score_list, names)))
score_list <- lapply(score_list, function(d) { for (cc in setdiff(all_cols, names(d))) d[[cc]] <- NA_real_; d[, all_cols, drop = FALSE] })
scores <- do.call(rbind, score_list)
gene_audit <- do.call(rbind, lapply(res, `[[`, "genes"))
write.csv(scores, file.path(out, "mechanism_module_scores.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(gene_audit, file.path(out, "mechanism_module_gene_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")

grp_summary <- do.call(rbind, lapply(split(scores, scores$dataset), function(d) {
  mods <- names(modules)
  do.call(rbind, lapply(mods, function(m) {
    gg <- as.character(d$group); keep <- !is.na(gg) & !is.na(d[[m]])
    if (!any(keep)) return(data.frame(dataset = d$dataset[1], module = m, group = NA_character_, median_score = NA_real_))
    ug <- sort(unique(gg[keep]))
    data.frame(dataset = d$dataset[1], module = m, group = ug,
               median_score = vapply(ug, function(g) median(d[[m]][keep & gg == g]), numeric(1)), stringsAsFactors = FALSE)
  }))
}))
write.csv(grp_summary, file.path(out, "mechanism_module_group_medians.csv"), row.names = FALSE, fileEncoding = "UTF-8")

pairwise <- do.call(rbind, lapply(split(scores, scores$dataset), function(d) {
  mods <- names(modules)
  do.call(rbind, lapply(mods, function(m) {
    if (d$dataset[1] == "GSE24206") { a <- d[[m]][d$group == "early_IPF"]; b <- d[[m]][d$group == "control"]; contrast <- "early_IPF_vs_control" }
    else { a <- d[[m]][d$group == "IPF_UIP" | d$group == "UIP"]; b <- d[[m]][d$group == "control"]; contrast <- "IPF_UIP_or_UIP_vs_control" }
    p <- if (length(a) > 1 && length(b) > 1) wilcox.test(a, b, exact = FALSE)$p.value else NA_real_
    data.frame(dataset = d$dataset[1], module = m, contrast = contrast, n_case = length(a), n_control = length(b), median_case = median(a), median_control = median(b), wilcox_p = p, stringsAsFactors = FALSE)
  }))
}))
pairwise$wilcox_fdr <- p.adjust(pairwise$wilcox_p, method = "BH")
write.csv(pairwise, file.path(out, "mechanism_module_pairwise_tests.csv"), row.names = FALSE, fileEncoding = "UTF-8")

capture.output({
  cat("Prespecified mechanism pilot; no reflux, LPR, aspiration, or laryngoscopy exposure is measured.\n\n")
  print(gene_audit)
  cat("\nGroup counts:\n"); print(with(scores, table(dataset, group)))
  cat("\nModule group medians:\n"); print(grp_summary)
  cat("\nExploratory pairwise tests (Wilcoxon; not reflux exposure tests):\n"); print(pairwise)
}, file = file.path(out, "mechanism_pilot_report.txt"))
