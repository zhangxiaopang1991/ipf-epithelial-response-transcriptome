const fs = require('fs');
const zlib = require('zlib');
const readline = require('readline');

const project = process.cwd();
const matrixPath = `${project}/inputs/GSE136831_RawCounts_Sparse.mtx.gz`;
const outPath = `${project}/outputs/GSE136831_cell_library_sizes.tsv`;
const dimsExpected = [45947, 312928, 692789348];
const totals = new Float64Array(dimsExpected[1]);
const input = readline.createInterface({ input: fs.createReadStream(matrixPath).pipe(zlib.createGunzip()), crlfDelay: Infinity });
let header = false, dims = null, lines = 0;
input.on('line', line => {
  if (!header) { if (line.startsWith('%%MatrixMarket')) header = true; return; }
  if (!dims) { if (line.startsWith('%') || !line.trim()) return; dims = line.trim().split(/\s+/).map(Number); if (dims.join(',') !== dimsExpected.join(',')) throw new Error(`Unexpected dimensions: ${dims}`); return; }
  const p = line.trim().split(/\s+/);
  if (p.length < 3) throw new Error(`Malformed matrix row at ${lines + 1}`);
  const col = Number(p[1]) - 1;
  const count = Number(p[2]);
  if (!Number.isInteger(col) || col < 0 || col >= dimsExpected[1] || !Number.isFinite(count) || count < 0) throw new Error(`Invalid entry at ${lines + 1}`);
  totals[col] += count;
  lines++;
  if (lines % 100000000 === 0) process.stderr.write(`scanned=${lines}\n`);
});
input.on('close', () => {
  if (lines !== dimsExpected[2]) throw new Error(`Entry count mismatch: ${lines}`);
  const out = fs.createWriteStream(outPath, { encoding: 'utf8' });
  out.write('cell_col\tlibrary_umi\n');
  for (let i = 0; i < totals.length; i++) out.write(`${i + 1}\t${totals[i]}\n`);
  out.end(() => process.stdout.write(JSON.stringify({ dims, scanned_entries: lines, cells: totals.length, output: outPath }) + '\n'));
});
input.on('error', err => { console.error(err.stack || err); process.exitCode = 1; });
