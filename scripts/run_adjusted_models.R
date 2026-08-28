library(GEOquery)

root <- normalizePath(".", mustWork = TRUE)
out <- normalizePath("outputs", mustWork = TRUE)
sc <- read.csv(file.path(out, "mechanism_module_scores.csv"), check.names = FALSE)
sc <- subset(sc, dataset == "GSE32537")
es <- getGEO(filename = file.path("inputs/geo_series_matrix", "GSE32537_series_matrix.txt.gz"), GSEMatrix = TRUE, getGPL = FALSE)
pd <- pData(es); rownames(pd) <- as.character(pd$geo_accession)
cc <- pd[, grep("^characteristics", names(pd)), drop = FALSE]
get_field <- function(pattern, numeric = FALSE) {
  z <- apply(cc, 1, function(x) { y <- x[grepl(pattern, tolower(x))]; if (length(y)) y[1] else NA_character_ })
  z <- sub("^[^:]*:\\s*", "", z)
  if (numeric) as.numeric(sub("[^0-9.-].*$", "", z)) else factor(z)
}
meta <- data.frame(sample_id = rownames(pd), age = get_field("^age", TRUE), gender = get_field("^gender"), smoking = get_field("smoking status"), rin = get_field("^rin", TRUE), tissue_source = get_field("tissue source"), stringsAsFactors = FALSE)
dat <- merge(sc, meta, by = "sample_id", all.x = TRUE, sort = FALSE)
dat$case <- factor(ifelse(dat$group == "IPF_UIP", "IPF_UIP", ifelse(dat$group == "control", "control", NA)), levels = c("control", "IPF_UIP"))
mods <- c("epithelial_barrier", "mucociliary_clearance", "injury_inflammation", "oxidative_stress", "fibrosis_epithelial_response")
fit_rows <- list(); i <- 0
for (m in mods) {
  d <- dat[!is.na(dat$case), c(m, "case", "age", "gender", "smoking", "rin", "tissue_source")]
  d <- d[complete.cases(d), , drop = FALSE]
  f <- as.formula(paste(m, "~ case + age + gender + smoking + rin + tissue_source"))
  fit <- lm(f, data = d); sm <- summary(fit)$coefficients
  term <- "caseIPF_UIP"
  i <- i + 1; fit_rows[[i]] <- data.frame(outcome = m, analysis = "IPF_UIP_vs_control_adjusted", n = nrow(d), beta = sm[term, 1], se = sm[term, 2], p = sm[term, 4], r2 = summary(fit)$r.squared, stringsAsFactors = FALSE)
}
for (v in c("FVC_percent", "DLCO_percent")) for (m in mods) {
  d <- dat[dat$group == "IPF_UIP", c(v, m, "age", "gender", "smoking", "rin", "tissue_source")]
  d <- d[complete.cases(d), , drop = FALSE]
  f <- as.formula(paste(v, "~", m, "+ age + gender + smoking + rin + tissue_source"))
  fit <- lm(f, data = d); sm <- summary(fit)$coefficients
  i <- i + 1; fit_rows[[i]] <- data.frame(outcome = v, analysis = paste0("IPF_UIP_module_association_", m), n = nrow(d), beta = sm[m, 1], se = sm[m, 2], p = sm[m, 4], r2 = summary(fit)$r.squared, stringsAsFactors = FALSE)
}
fits <- do.call(rbind, fit_rows); fits$fdr_within_all_models <- p.adjust(fits$p, method = "BH")
write.csv(fits, file.path(out, "adjusted_linear_models_GSE32537.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# Sensitivity analysis omitting tissue source, whose levels are strongly collinear with other covariates.
sens_rows <- list(); j <- 0
for (m in mods) {
  d <- dat[!is.na(dat$case), c(m, "case", "age", "gender", "smoking", "rin")]; d <- d[complete.cases(d), , drop = FALSE]
  fit <- lm(as.formula(paste(m, "~ case + age + gender + smoking + rin")), d); sm <- summary(fit)$coefficients; j <- j + 1
  sens_rows[[j]] <- data.frame(outcome = m, analysis = "IPF_UIP_vs_control_sensitivity_no_tissue_source", n = nrow(d), beta = sm["caseIPF_UIP",1], se = sm["caseIPF_UIP",2], p = sm["caseIPF_UIP",4], r2 = summary(fit)$r.squared)
}
sens <- do.call(rbind, sens_rows); sens$fdr <- p.adjust(sens$p, "BH")
write.csv(sens, file.path(out, "adjusted_model_sensitivity_no_tissue_source.csv"), row.names = FALSE, fileEncoding = "UTF-8")
capture.output({ cat("GSE32537 adjusted exploratory linear models; no reflux exposure measured.\n"); print(fits) }, file = file.path(out, "adjusted_models_report.txt"))
