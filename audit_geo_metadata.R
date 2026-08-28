base <- normalizePath(file.path("inputs", "geo_series_matrix"), mustWork = TRUE)
out <- normalizePath("outputs", mustWork = TRUE)

parse_quoted <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) < 2) return(character())
  gsub('"', '', fields[-1], fixed = TRUE)
}
parse_fields <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(fields) < 2) return(character())
  gsub('"', '', fields[-1], fixed = TRUE)
}

files <- list.files(base, pattern = "_series_matrix\\.txt\\.gz$", full.names = TRUE)
rows <- list()
for (f in files) {
  con <- gzfile(f, "rt")
  on.exit(close(con), add = TRUE)
  sample_ids <- character(); titles <- character(); sources <- character(); chars <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (!length(line)) break
    if (startsWith(line, "!Series_sample_id")) {
      id_field <- parse_quoted(line)
      sample_ids <- unlist(strsplit(trimws(gsub('"', '', id_field[1], fixed = TRUE)), " +"))
    }
    if (startsWith(line, "!Sample_title")) titles <- parse_fields(line)
    if (startsWith(line, "!Sample_source_name_ch1")) sources <- parse_fields(line)
    if (startsWith(line, "!Sample_characteristics_ch1")) chars[[length(chars) + 1]] <- parse_fields(line)
    if (startsWith(line, "!series_matrix_table_begin")) break
  }
  close(con)
  gse <- sub("_series_matrix.*$", "", basename(f))
  n <- length(sample_ids)
  if (!n) next
  if (length(titles) && length(titles) != n) stop(sprintf("%s: title count %d != sample count %d", gse, length(titles), n))
  if (length(sources) && length(sources) != n) stop(sprintf("%s: source count %d != sample count %d", gse, length(sources), n))
  if (length(chars) && any(vapply(chars, length, integer(1)) != n)) stop(sprintf("%s: characteristic field length mismatch", gse))
  char_text <- if (length(chars)) vapply(seq_len(n), function(i) paste(vapply(chars, `[`, "", i), collapse = " | "), "") else rep("", n)
  rows[[length(rows) + 1]] <- data.frame(
    dataset = gse, sample_id = sample_ids,
    title = if(length(titles)==n) titles else rep(NA_character_, n),
    source = if(length(sources)==n) sources else rep(NA_character_, n),
    characteristics = char_text, stringsAsFactors = FALSE
  )
}
audit <- do.call(rbind, rows)
write.csv(audit, file.path(out, "geo_sample_metadata_audit.csv"), row.names = FALSE, fileEncoding = "UTF-8")
summary <- do.call(rbind, lapply(split(audit, audit$dataset), function(d) {
  data.frame(dataset = d$dataset[1], samples = nrow(d),
             ipf_terms = sum(grepl("IPF|UIP|fibrosis", paste(d$title, d$source, d$characteristics), ignore.case = TRUE)),
             reflux_terms = sum(grepl("reflux|aspirat|GERD|LPR", paste(d$title, d$source, d$characteristics), ignore.case = TRUE)),
             fvc_terms = sum(grepl("FVC|forced vital capacity", d$characteristics, ignore.case = TRUE)),
             dlco_terms = sum(grepl("DLCO|DLCO", d$characteristics, ignore.case = TRUE)),
             stringsAsFactors = FALSE)
}))
write.csv(summary, file.path(out, "geo_dataset_audit_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
capture.output({ print(summary); cat("\\nNo reflux exposure field is inferred from absence of matching terms.\\n") }, file = file.path(out, "geo_dataset_audit_summary.txt"))
