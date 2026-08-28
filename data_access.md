# Data access and provenance

The analysis uses public NCBI GEO accessions GSE10667, GSE24206, and GSE32537, with GPL4133, GPL570, and GPL6244 platform annotations. The audited input checksums are in `data_provenance.md`.

## Manual download and placement

Download the series-matrix archive from each GEO accession page, then place the unmodified files at:

- `inputs/geo_series_matrix/GSE10667_series_matrix.txt.gz`
- `inputs/geo_series_matrix/GSE24206_series_matrix.txt.gz`
- `inputs/geo_series_matrix/GSE32537_series_matrix.txt.gz`

Download each official GPL annotation archive and place the unmodified files at:

- `inputs/GPL4133.annot.gz`
- `inputs/GPL570.annot.gz`
- `inputs/GPL6244.annot.gz`

The `inputs/` directory is excluded from Git. Do not silently alter source files. Public GEO sample labels may be retained only where required for reproducibility; do not infer or add donor identities.
