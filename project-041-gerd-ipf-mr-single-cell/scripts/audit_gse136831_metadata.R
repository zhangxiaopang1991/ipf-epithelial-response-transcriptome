options(stringsAsFactors = FALSE)

project <- normalizePath(".")
soft <- file.path(project, "inputs", "GSE136831_family.soft.gz")
meta_file <- file.path(project, "inputs", "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz")
out <- file.path(project, "outputs")

z <- gzfile(soft, "rt")
x <- readLines(z, warn = FALSE)
close(z)
starts <- which(grepl("^\\^SAMPLE", x))
ends <- c(starts[-1] - 1L, length(x))
sample_rows <- lapply(seq_along(starts), function(i) {
  b <- x[starts[i]:ends[i]]
  field <- function(pattern) {
    v <- b[grepl(pattern, b)]
    if (!length(v)) return(NA_character_)
    sub(pattern, "", v[1])
  }
  data.frame(
    geo_accession = field("^\\^SAMPLE = "),
    title = field("^!Sample_title = "),
    disease = field("^!Sample_characteristics_ch1 = disease: "),
    sex = field("^!Sample_characteristics_ch1 = Sex: "),
    age = field("^!Sample_characteristics_ch1 = age: "),
    source_name = field("^!Sample_source_name_ch1 = "),
    stringsAsFactors = FALSE
  )
})
series_samples <- do.call(rbind, sample_rows)
write.table(series_samples, file.path(out, "GSE136831_series_sample_metadata.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, na = "")

con <- gzfile(meta_file, "rt")
cell <- read.delim(con, check.names = FALSE, quote = "\"", stringsAsFactors = FALSE)
close(con)
names(cell) <- gsub('"', '', names(cell), fixed = TRUE)
for (j in seq_along(cell)) cell[[j]] <- gsub('"', '', cell[[j]], fixed = TRUE)

series_samples$subject_from_title <- sub(" scRNAseq$", "", series_samples$title)
subject_map <- unique(cell[c("Subject_Identity", "Disease_Identity")])
mapping_check <- merge(series_samples[c("geo_accession", "title", "disease", "subject_from_title")],
                       subject_map, by.x = "subject_from_title", by.y = "Subject_Identity", all.x = TRUE)
mapping_check$subject_present_in_cell_metadata <- !is.na(mapping_check$Disease_Identity)
mapping_check$disease_matches <- mapping_check$disease == mapping_check$Disease_Identity
write.table(mapping_check, file.path(out, "GSE136831_series_subject_mapping_check.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

group_summary <- do.call(rbind, lapply(split(cell, cell$Disease_Identity), function(d) {
  data.frame(
    disease = unique(d$Disease_Identity),
    cells = nrow(d),
    subjects = length(unique(d$Subject_Identity)),
    libraries = length(unique(d$Library_Identity)),
    cell_type_category_n = length(unique(d$CellType_Category)),
    manuscript_identity_n = length(unique(d$Manuscript_Identity)),
    subclass_identity_n = length(unique(d$Subclass_Cell_Identity)),
    stringsAsFactors = FALSE
  )
}))
write.table(group_summary, file.path(out, "GSE136831_cell_metadata_group_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

subject_summary <- aggregate(CellBarcode_Identity ~ Disease_Identity + Subject_Identity, cell, length)
names(subject_summary)[3] <- "cells"
subject_summary$libraries <- vapply(seq_len(nrow(subject_summary)), function(i) {
  d <- cell[cell$Disease_Identity == subject_summary$Disease_Identity[i] & cell$Subject_Identity == subject_summary$Subject_Identity[i], ]
  length(unique(d$Library_Identity))
}, integer(1))
write.table(subject_summary, file.path(out, "GSE136831_subject_cell_counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

celltype_summary <- as.data.frame(table(cell$Disease_Identity, cell$Subclass_Cell_Identity), stringsAsFactors = FALSE)
names(celltype_summary) <- c("disease", "subclass_cell_identity", "cells")
celltype_summary <- celltype_summary[celltype_summary$cells > 0, ]
write.table(celltype_summary, file.path(out, "GSE136831_celltype_by_disease.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

cat("Series samples:", nrow(series_samples), "\\n")
cat("Metadata cells:", nrow(cell), "\\n")
cat("Metadata subjects by disease:\\n")
print(with(cell, table(Disease_Identity, Subject_Identity)))
cat("Metadata libraries by disease:\\n")
print(with(cell, tapply(Library_Identity, Disease_Identity, function(v) length(unique(v)))))
cat("Missingness:\\n")
print(colSums(is.na(cell) | cell == ""))
