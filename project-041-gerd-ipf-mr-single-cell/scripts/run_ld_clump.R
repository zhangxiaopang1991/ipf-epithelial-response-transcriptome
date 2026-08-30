args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: run_ld_clump.R input.tsv output.tsv")
.libPaths(c(file.path(normalizePath("."), ".r_libs"), .libPaths()))
script_arg <- commandArgs()[grep("^--file=", commandArgs())]
if (length(script_arg)) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/")
  project_root <- dirname(dirname(dirname(dirname(script_path))))
  .libPaths(c(file.path(project_root, ".r_libs"), .libPaths()))
}
library(ieugwasr)

input <- args[[1]]
output <- args[[2]]
d <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE)
rs <- if ("rsid" %in% names(d)) d$rsid else d$hm_rsid
dat <- data.frame(rsid = rs, pval = d$p_value, id = "trait")
dat <- dat[!is.na(dat$rsid) & dat$rsid != "", ]
jwt <- Sys.getenv("OPENGWAS_JWT", unset = "")
if (!nzchar(jwt)) stop("OPENGWAS_JWT is not set; no token is written to disk")
clumped <- ieugwasr::ld_clump(dat, clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8, pop = "EUR", opengwas_jwt = jwt)
write.table(clumped, output, sep = "\t", row.names = FALSE, quote = FALSE)
cat("retained", nrow(clumped), "variants\n")
