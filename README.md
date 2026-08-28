# Public code package draft (v1)

This folder is a release candidate for a public repository accompanying the IPF lung-tissue transcriptome manuscript. It intentionally contains no GEO expression matrices, platform annotation archives, participant-level metadata, author information, or credentials.

## Intended repository contents

- `scripts/`: the analysis R scripts copied from the audited project `work/` directory.
- `run_order.md`: the documented execution order and expected outputs.
- `software_versions.md`: R/Bioconductor package requirements.
- `data_access.md`: public GEO accessions and download/checksum policy.
- `disclosure_boundary.md`: what is and is not released.

## Before public release

1. Copy this folder to a repository selected by the authors.
2. Add a tested download script or precise manual download instructions for GSE10667, GSE24206, GSE32537 and GPL4133/GPL570/GPL6244.
3. Add a clean `sessionInfo()` output generated from the release environment.
4. Run the complete order from a clean checkout and compare declared output assertions.
5. Review whether GEO terms-of-use and the target journal require a persistent archive DOI.

This draft is not itself a public repository and must not be cited as one until released.
