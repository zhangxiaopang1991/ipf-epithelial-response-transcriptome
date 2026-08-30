# Software and dependency requirements

The analysis was developed with R 4.6.1, Node.js 20.11.0 and the following R packages: data.table 1.18.4, MendelianRandomization 0.10.0, MRPRESSO 1.0, ieugwasr 1.1.0.9000, ggplot2 4.0.3 and R.utils 2.13.0.

Install the required R packages in a user library before running the scripts. The OpenGWAS clumping step additionally requires a valid `OPENGWAS_JWT` environment variable; the token must never be written to files or committed.

The repository intentionally excludes large GWAS and single-cell input files. Download them from the accessions documented in `README.md`, verify their checksums where provided, and place them in the project-relative input locations expected by the scripts.

The exact package versions above describe the analysis environment. Minor-version differences may alter numerical output and should be recorded in the reproduction log.
