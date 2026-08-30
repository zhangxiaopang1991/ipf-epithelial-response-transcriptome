$ErrorActionPreference = 'Stop'
$project = (Get-Location).Path
$inputFile = Join-Path $project 'outputs\GERD_GCST90000514_ld_clumped.tsv'
$outputFile = Join-Path $project 'outputs\GERD_GCST90000514_ld_clumped_ensembl_vep.tsv'
$rsids = Import-Csv $inputFile -Delimiter "`t" | ForEach-Object { $_.rsid }
$body = @{ ids = @($rsids) } | ConvertTo-Json -Depth 4
$result = Invoke-RestMethod -Uri 'https://rest.ensembl.org/vep/human/id' -Method Post -Body $body -ContentType 'application/json' -Headers @{ Accept = 'application/json' }
$rows = foreach ($v in $result) {
  $id = if ($v.id) { $v.id } else { $v.input }
  $genes = @($v.transcript_consequences | Where-Object { $_.gene_symbol } | ForEach-Object { $_.gene_symbol } | Sort-Object -Unique)
  [pscustomobject]@{
    rsid = $id
    assembly = $v.assembly_name
    chromosome = $v.seq_region_name
    position = $v.start
    most_severe_consequence = $v.most_severe_consequence
    gene_symbols = ($genes -join ';')
    gene_count = $genes.Count
  }
}
$rows = $rows | Group-Object rsid | ForEach-Object {
  $g = $_.Group
  $genes = @($g.gene_symbols -split ';' | Where-Object { $_ } | Sort-Object -Unique)
  [pscustomobject]@{
    rsid = $_.Name
    assembly = ($g.assembly | Where-Object { $_ } | Select-Object -First 1)
    chromosome = ($g.chromosome | Where-Object { $_ } | Select-Object -First 1)
    position = ($g.position | Where-Object { $_ } | Select-Object -First 1)
    most_severe_consequence = (($g.most_severe_consequence | Sort-Object -Unique) -join ';')
    gene_symbols = ($genes -join ';')
    gene_count = $genes.Count
  }
}
$rows | Sort-Object rsid | Export-Csv $outputFile -Delimiter "`t" -NoTypeInformation -Encoding utf8
Write-Output ("annotated={0}; output={1}" -f $rows.Count, $outputFile)
