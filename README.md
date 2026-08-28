# Public code package (v1.0.0)

This folder contains the reproducible analysis code accompanying the IPF lung-tissue transcriptome manuscript. It intentionally contains no GEO expression matrices, platform annotation archives, participant-level metadata beyond public GEO sample labels, author information, or credentials.

## Intended repository contents

- `scripts/`: the analysis R scripts copied from the audited project `work/` directory.
- `run_order.md`: the documented execution order and expected outputs.
- `software_versions.md`: R/Bioconductor package requirements.
- `data_access.md`: public GEO accessions and download/checksum policy.
- `data_provenance.md`: source URLs, access date, and SHA-256 checksums for the inputs used for the audited run.
- `disclosure_boundary.md`: what is and is not released.

## Reproduction

1. Download the public GEO series-matrix and GPL annotation files using the exact instructions in `data_access.md`.
2. Place them in the exact `inputs/` paths described there. The directory is ignored by Git.
3. Run the documented commands in `run_order.md` from the repository root.
4. Compare generated outputs against the acceptance criteria in `reproduction_acceptance_criteria.md`.

The reported analysis passed a clean-directory reproduction audit on 2026-08-28. Hallmark GSEA is optional and is not reported in the manuscript.

## Release scope

This repository is code and documentation only. It does not make inferences about reflux, laryngopharyngeal reflux, microaspiration, laryngoscopy, causality, prognosis, or prediction.
