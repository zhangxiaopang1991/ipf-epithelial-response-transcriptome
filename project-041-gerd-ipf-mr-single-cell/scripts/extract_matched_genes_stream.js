const fs = require('fs');
const zlib = require('zlib');
const readline = require('readline');
const path = require('path');

const project = process.cwd();
const geneIndexPath = path.join(project, 'inputs/GSE136831_AllCells.GeneIDs.txt.gz');
const matchPath = path.join(project, 'outputs/GERD_candidate_genes_GSE136831_index_match.tsv');
const matrixPath = path.join(project, 'inputs/GSE136831_RawCounts_Sparse.mtx.gz');
const outPath = path.join(project, 'outputs/GSE136831_matched_genes_sparse.tsv');

const geneText = zlib.gunzipSync(fs.readFileSync(geneIndexPath)).toString('utf8').trim().split(/\r?\n/);
const geneRows = geneText.slice(1).map(x => x.split('\t').map(v => v.replace(/^"|"$/g, '')));
const wanted = new Map();
const matchLines = fs.readFileSync(matchPath, 'utf8').trim().split(/\r?\n/).slice(1);
for (const line of matchLines) {
  const [symbol, ensembl, matched] = line.split('\t');
  if (matched === 'TRUE' && ensembl) {
    const row = geneRows.findIndex(x => x[0] === ensembl) + 1;
    if (row > 0) wanted.set(row, { ensembl, symbol });
  }
}

const input = readline.createInterface({ input: fs.createReadStream(matrixPath).pipe(zlib.createGunzip()), crlfDelay: Infinity });
const out = fs.createWriteStream(outPath, { encoding: 'utf8' });
out.write('gene_row\tcell_col\tcount\tgene_id\tgene_symbol\n');
let headerSeen = false;
let dims = null;
let kept = 0;
let lines = 0;
input.on('line', line => {
  if (!headerSeen) {
    if (line.startsWith('%%MatrixMarket')) headerSeen = true;
    return;
  }
  if (!dims) {
    if (line.startsWith('%') || !line.trim()) return;
    dims = line.trim().split(/\s+/).map(Number);
    return;
  }
  lines++;
  const p = line.trim().split(' ');
  const row = Number(p[0]);
  const hit = wanted.get(row);
  if (hit) {
    const col = Number(p[1]);
    const count = Number(p[2]);
    out.write(`${row}\t${col}\t${count}\t${hit.ensembl}\t${hit.symbol}\n`);
    kept++;
  }
  if (lines % 100000000 === 0) process.stderr.write(`scanned=${lines}; kept=${kept}\n`);
});
input.on('close', () => {
  out.end(() => {
    process.stdout.write(JSON.stringify({ dims, target_genes: wanted.size, scanned_entries: lines, kept_entries: kept, output: outPath }) + '\n');
  });
});
