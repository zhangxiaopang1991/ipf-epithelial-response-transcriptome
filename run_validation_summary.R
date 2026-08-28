library(GEOquery)

out <- normalizePath("outputs", mustWork = TRUE)
fits <- read.csv(file.path(out, "adjusted_linear_models_GSE32537.csv"))
sc <- read.csv(file.path(out, "mechanism_module_scores.csv"), check.names = FALSE)
pd <- pData(getGEO(filename = file.path("inputs/geo_series_matrix", "GSE32537_series_matrix.txt.gz"), GSEMatrix = TRUE, getGPL = FALSE))
cc <- pd[, grep("^characteristics", names(pd)), drop = FALSE]
get_field <- function(pattern, numeric = FALSE) { z <- apply(cc, 1, function(x) { y <- x[grepl(pattern, tolower(x))]; if (length(y)) y[1] else NA_character_ }); z <- sub("^[^:]*:\\s*", "", z); if (numeric) as.numeric(sub("[^0-9.-].*$", "", z)) else factor(z) }
meta <- data.frame(sample_id = pd$geo_accession, age = get_field("^age", TRUE), gender = get_field("^gender"), smoking = get_field("smoking status"), rin = get_field("^rin", TRUE), tissue_source = get_field("tissue source"))
dat <- merge(subset(sc, dataset == "GSE32537"), meta, by = "sample_id", sort = FALSE)
mods <- c("epithelial_barrier", "mucociliary_clearance", "injury_inflammation", "oxidative_stress", "fibrosis_epithelial_response")
diag <- do.call(rbind, lapply(mods, function(m) { d <- dat[!is.na(dat$group) & dat$group %in% c("IPF_UIP", "control"), ]; d <- d[complete.cases(d[, c(m,"group","age","gender","smoking","rin","tissue_source")]), ]; fit <- lm(as.formula(paste(m,"~ group + age + gender + smoking + rin + tissue_source")), d); mm <- model.matrix(fit); v <- vapply(seq_len(ncol(mm)-1)+1, function(j) { r2 <- summary(lm(mm[,j] ~ mm[,-j]))$r.squared; 1/(1-r2) }, numeric(1)); h <- hatvalues(fit); stdres <- residuals(fit)/(summary(fit)$sigma*sqrt(pmax(1-h, .Machine$double.eps))); cook <- (stdres^2/length(coef(fit)))*(h/pmax(1-h,.Machine$double.eps)); data.frame(outcome=m, n=nrow(d), residual_sd=sd(residuals(fit)), max_abs_std_resid=max(abs(stdres),na.rm=TRUE), max_cooks=max(cook,na.rm=TRUE), max_vif=max(v), stringsAsFactors=FALSE) }))
write.csv(diag, file.path(out, "adjusted_model_diagnostics.csv"), row.names = FALSE, fileEncoding = "UTF-8")

p <- read.csv(file.path(out, "mechanism_module_pairwise_tests.csv"))
p$delta <- p$median_case - p$median_control
p$direction <- ifelse(p$delta > 0, "higher_in_case", ifelse(p$delta < 0, "lower_in_case", "zero"))
conc <- do.call(rbind, lapply(split(p, p$module), function(d) data.frame(module=d$module[1], datasets=nrow(d), positive_delta=sum(d$delta>0), negative_delta=sum(d$delta<0), direction_pattern=paste(paste(d$dataset,d$direction,sep=":"),collapse=";"), stringsAsFactors=FALSE)))
write.csv(conc, file.path(out, "cross_dataset_direction_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
capture.output({ cat("Adjusted model diagnostics\n"); print(diag); cat("\nCross-dataset direction summary\n"); print(conc) }, file = file.path(out, "validation_summary_report.txt"))
