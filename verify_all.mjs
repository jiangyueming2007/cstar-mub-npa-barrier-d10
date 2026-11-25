#!/usr/bin/env node
// One-click ZERO-FLOAT verification suite for all core theorems (Node.js v20+, ZERO external dependencies!).
// Usage:
//   node verify_all.mjs          -> runs all 8 zero-float Q.E.D. proofs (~2 minutes total)
//   node verify_all.mjs --paper2 -> runs only Paper II (C*-MUB / NPA Barrier) suite (~1.5 minutes)
//   node verify_all.mjs --fast   -> runs fast core proofs (<35s total)

import { execFileSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { fileURLToPath } from 'node:url';

const __dir = path.dirname(fileURLToPath(import.meta.url));
const mode = process.argv[2] || '--all';

const SUITE = [
  { id: '01', paper: 2, fast: true,  file: '01_quaternionic_4_5_mub_QED.mjs',              desc: 'Quaternionic 4-MUB (d=6) & 5-MUB (d=10,14) over Z[i,w_p][j]' },
  { id: '02', paper: 2, fast: true,  file: '02_complex_M12_M20_projectors_QED.mjs',        desc: 'Explicit M_12 & M_20 rank-2 projectors (432 + 2,000 triple checks)' },
  { id: '03', paper: 2, fast: true,  file: '03_gribling_polak_def31_alignment_QED.mjs',    desc: 'Gribling-Polak Def. 3.1 axioms (i)-(iv) + deg-7/10 witnesses' },
  { id: '04', paper: 2, fast: false, file: '04_clifford_6_8_12_mub_QED.mjs',               desc: 'Clifford (q+1)-MUBs: d=10 6-MUB (M_40), d=14 8-MUB, d=22 12-MUB' },
  { id: '05', paper: 2, fast: true,  file: '05_universal_jordan_wigner_N_Cstar_inf_QED.mjs', desc: 'Universal Jordan-Wigner Theorem: N_{C*}(d) = infinity for all d >= 2' },
  { id: '06', paper: 1, fast: true,  file: '06_B_new_no_4th_vec_nullstellensatz_QED.mjs',  desc: 'B_new (36 bases) Degree-6 Nullstellensatz (1256/1256 mod p1, p2)' },
  { id: '07', paper: 1, fast: true,  file: '07_B_new_fixed_Hconf_isolation_QED.mjs',       desc: 'B_new fixed-{I, H_conf} 90/90 isolation & M_4 4/4 transversality' },
  { id: '08', paper: 1, fast: true,  file: '08_B_iso_72orbits_171D_rigidity_QED.mjs',      desc: 'B_iso 72 chiral orbits & 171/171 Joint 3-MUB Rigidity' },
];

const selected = SUITE.filter(s => {
  if (mode === '--paper2') return s.paper === 2;
  if (mode === '--paper1') return s.paper === 1;
  if (mode === '--fast') return s.fast;
  return true;
});

console.log(`=== ZERO-FLOAT MUB VERIFICATION SUITE (mode=${mode}, ${selected.length} proofs) ===\n`);
const tAll = Date.now();
const summary = [];

for (const item of selected) {
  const fullPath = path.join(__dir, 'scripts', item.file);
  const sha = crypto.createHash('sha256').update(fs.readFileSync(fullPath)).digest('hex').slice(0, 12);
  console.log(`[${item.id}] Running ${item.file} (sha256:${sha}...) — ${item.desc}`);
  const t0 = Date.now();
  execFileSync(process.execPath, [fullPath], { stdio: 'inherit' });
  const dt = +((Date.now() - t0) / 1000).toFixed(2);
  summary.push({ id: item.id, file: item.file, sha256_prefix: sha, sec: dt, status: 'Q.E.D. PASS' });
  console.log(`    --> PASS in ${dt}s\n`);
}

console.log('=== ALL ZERO-FLOAT PROOFS PASSED ===');
for (const s of summary) console.log(`  [${s.id}] ${s.status} | ${String(s.sec).padStart(6)}s | ${s.sha256_prefix} | ${s.file}`);
console.log(`Total wall time: ${((Date.now() - tAll) / 1000).toFixed(2)}s`);
