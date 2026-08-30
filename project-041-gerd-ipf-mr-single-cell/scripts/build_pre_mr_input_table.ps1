$ErrorActionPreference = 'Stop'
$project = (Get-Location).Path
$harm = Import-Csv (Join-Path $project 'outputs\GERD_to_IPF_harmonisation_check.tsv') -Delimiter "`t"
$ann = Import-Csv (Join-Path $project 'outputs\GERD_GCST90000514_ld_clumped_ensembl_vep.tsv') -Delimiter "`t"
$annById = @{}
foreach ($a in $ann) { $annById[$a.rsid] = $a }
$rows = foreach ($h in $harm) {
  $f = ([double]$h.exposure_beta / [double]$h.exposure_se) * ([double]$h.exposure_beta / [double]$h.exposure_se)
  $a = $annById[$h.rsid]
  [pscustomobject]@{
    rsid = $h.rsid; exposure_p = $h.exposure_p; exposure_beta = $h.exposure_beta; exposure_se = $h.exposure_se; exposure_eaf = $h.exposure_eaf
    outcome_beta = $h.outcome_beta; outcome_se = $h.outcome_se; outcome_eaf = $h.outcome_eaf; palindromic = $h.palindromic; palindromic_ambiguous = $h.palindromic_ambiguous
    allele_status = $h.allele_status; F_statistic = [math]::Round($f, 4); primary_eligible = ($h.palindromic_ambiguous -ne 'True'); gene_symbols = $a.gene_symbols; most_severe_consequence = $a.most_severe_consequence
  }
}
$out = Join-Path $project 'outputs\GERD_to_IPF_pre_mr_input_table.tsv'
$rows | Sort-Object rsid | Export-Csv $out -Delimiter "`t" -NoTypeInformation -Encoding utf8
$fvals = @($rows.F_statistic | Sort-Object)
$mid = [math]::Floor($fvals.Count / 2)
$median = ($fvals[$mid - 1] + $fvals[$mid]) / 2
[pscustomobject]@{ total = $rows.Count; primary_eligible = @($rows | Where-Object primary_eligible).Count; ambiguous_palindromic_excluded = @($rows | Where-Object { -not $_.primary_eligible }).Count; median_F = [math]::Round($median, 4); min_F = ($fvals | Measure-Object -Minimum).Minimum } | ConvertTo-Json
