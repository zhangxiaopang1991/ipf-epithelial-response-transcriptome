$project = (Get-Location).Path
$ann = Import-Csv (Join-Path $project 'outputs\GERD_GCST90000514_ld_clumped_ensembl_vep.tsv') -Delimiter "`t"
$genes = @($ann.gene_symbols -split ';' | Where-Object { $_ } | Sort-Object -Unique)
$out = Join-Path $project 'outputs\GERD_direct_annotated_candidate_genes.txt'
$genes | Set-Content $out -Encoding utf8
[pscustomobject]@{ snp_rows = $ann.Count; direct_gene_count = $genes.Count; output = $out } | ConvertTo-Json
