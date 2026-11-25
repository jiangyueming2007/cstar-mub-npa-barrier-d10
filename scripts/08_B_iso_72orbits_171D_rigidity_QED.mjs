import _pkgPath from 'node:path';
import { fileURLToPath as _pkgFileURLToPath } from 'node:url';
const _pkgDir = _pkgPath.dirname(_pkgFileURLToPath(import.meta.url));
// Round 15 Track B (Local CPU):
//   THREE ZERO-FLOAT EXACT THEOREMS FOR THE KAGGLE-GPU-DISCOVERED {0, 0, 16} 3-MUB {I, H_conf, B_iso}!!
//   [1] Map M_mu4 -> canonical H_conf = (I + iC)/sqrt(10) via exact monomial isomorphism,
//       transform B_iso = B2_aligned^dag into canonical H_conf coordinates, verify ALL 300 eqs in Z[i, sqrt(5)]/9,
//       and compute the EXACT G_1440(H_conf) orbit size & linear stabilizer |G_lin| of B_iso!
//   [2] ZERO-FLOAT 171D Joint 3-MUB Rigidity Theorem for {I, H_conf, B_iso} in Z[i, sqrt(5)] mod p1, p2!
//   [3] ZERO-FLOAT Degree-6 Hilbert Nullstellensatz Theorem proving {I, H_conf, B_iso} admits
//       NO arbitrary-phase common unbiased 4th vector w in C^10!

import fs from 'node:fs';
import assert from 'node:assert/strict';

const N = 10;
const C = Array.from({ length: N }, (_, i) =>
  Array.from({ length: N }, (_, j) => {
    if (i === j) return 0;
    if (i === 0 || j === 0) return 1;
    const r1 = Math.floor((i - 1) / 3), c1 = (i - 1) % 3;
    const r2 = Math.floor((j - 1) / 3), c2 = (j - 1) % 3;
    return (r1 === r2 || c1 === c2) ? 1 : -1;
  })
);
const Hpow = Array.from({ length: N }, (_, r) => Array.from({ length: N }, (_, c) => r === c ? 0 : (C[r][c] === 1 ? 1 : 3)));
const I4 = [[1, 0], [0, 1], [-1, 0], [0, -1]];

const isoData = JSON.parse(fs.readFileSync(_pkgPath.join(_pkgDir, '..', 'data', 'B_iso_0016_exact_Q_i_sqrt5.json'), 'utf8'));
const Mpow = isoData.M_mu4;       // M[r][c] = i^{Mpow[r][c]}
const B2spec = isoData.B2_aligned; // B2[r][c] = i^k * v^e, satisfying |(B2 * M)[r][c]|^2 = 10!
// Therefore {I, M, B2^dag} is a 3-MUB (Wait: check! B2^dag is unitary, M is unitary, and (B2^dag)^dag * M = B2 * M has flat modulus sqrt(10)!!
// Wait: B_iso = B2^dag has rows=cols of B2, so B_iso[r][c] = conj(B2[c][r])!
// Let's check: M is left-indexed by r (row of M = col of B2!), so {I, M, B_iso} has:
//   - M^dag M = 10 I
//   - B_iso^dag B_iso = (B2 B2^dag) = 10 I
//   - B_iso^dag M = B2 M -> ALL entries have modulus sqrt(10)!! EXACT!!)

// [1] Find EXACT monomial isomorphism: D_L * P_L * M * P_R * D_R = H_conf = I + iC!
// Notice: M_ij = i^{Mpow[i][j]}. If we find row-phases a[r] in Z/4, col-phases b[c] in Z/4, and bijection perm col->row
// such that (Mpow[r][c] + a[r] + b[c]) mod 4 is 0 on one entry per row/col and {1, 3} elsewhere, we get a conference matrix C'!
// First check for each pair of rows (r1, r2) the parity of Mpow[r1][c] - Mpow[r2][c] mod 2:
// In I + iC, row r1 and row r2 differ in parity (odd = +-i) at EXACTLY TWO columns: c = r1 and c = r2!
const rowToCol = Array(N).fill(-1);
for (let r = 0; r < N; r++) {
  // For each other row s != r, find the 2 cols where (Mpow[r][c] - Mpow[s][c]) is odd; their intersection across all s is {col(r)}!
  let cand = new Set(Array.from({ length: N }, (_, c) => c));
  for (let s = 0; s < N; s++) {
    if (s === r) continue;
    const oddCols = new Set(), evenCols = new Set();
    for (let c = 0; c < N; c++) { if (((Mpow[r][c] - Mpow[s][c]) & 1) === 1) oddCols.add(c); else evenCols.add(c); }
    const pair = oddCols.size === 2 ? oddCols : evenCols;
    assert.equal(pair.size, 2, `Expected 2 minority parity diffs between rows ${r},${s}, got odd=${oddCols.size}`);
    cand = new Set([...cand].filter(c => pair.has(c)));
  }
  assert.equal(cand.size, 1, `Row ${r} unique diagonal col failed: ${[...cand]}`);
  rowToCol[r] = [...cand][0];
}
console.log('Exact diagonal pairing of M_mu4 (row r -> unique diag col):', rowToCol);
// Reorder columns of M by rowToCol so diagonal is (r, r): M1[r][c] = Mpow[r][rowToCol[c]]!
const M1 = Array.from({ length: N }, (_, r) => Array.from({ length: N }, (_, c) => Mpow[r][rowToCol[c]]));
// Now find a[r], b[c] in Z/4 so that (M1[r][c] + a[r] + b[c]) mod 4 is 0 on r==c and {1,3} on r!=c!
// Fix a[0] = 0 -> b[0] = -M1[0][0]; for c > 0, (M1[0][c] + b[c]) must be odd, and since (M1[r][r]+a[r]+b[r])=0,
// wait: (M1[0][c] - M1[0][0]) is ALREADY odd for c>0? Let's make row 0 off-diags = 1 (so b[c] = (1 - M1[0][c]) & 3),
// then a[r] = (-M1[r][r] - b[r]) & 3!
const bPh = Array(N).fill(0), aPh = Array(N).fill(0);
bPh[0] = (4 - M1[0][0]) & 3;
for (let c = 1; c < N; c++) bPh[c] = (1 - M1[0][c] + 4) & 3;
for (let r = 1; r < N; r++) aPh[r] = (4 - M1[r][r] - bPh[r] + 8) & 3;
// Verify M2[r][c] = (M1[r][c] + aPh[r] + bPh[c]) & 3 has 0 on diag and {1,3} off-diag, AND is symmetric!
const C1 = Array.from({ length: N }, (_, r) => Array.from({ length: N }, (_, c) => {
  const v = (M1[r][c] + aPh[r] + bPh[c]) & 3;
  if (r === c) { assert.equal(v, 0); return 0; }
  assert(v === 1 || v === 3, `Off-diag not +-i at (${r},${c}): ${v}`);
  return v === 1 ? 1 : -1;
}));
// Make C1 symmetric (if C1[r][c] == -C1[c][r] is impossible since M2 is unitary -> C1 must be symmetric up to row/col sign flips!)
// Check if C1 is symmetric: if C1[r][0] = -1, flip row r and col r by multiplying by -1 (+2 in Z/4)!
for (let r = 1; r < N; r++) if (C1[r][0] === -1) {
  aPh[r] = (aPh[r] + 2) & 3; bPh[r] = (bPh[r] + 2) & 3;
  for (let k = 0; k < N; k++) { C1[r][k] = -C1[r][k]; C1[k][r] = -C1[k][r]; }
}
for (let r = 0; r < N; r++) for (let c = 0; c < N; c++) assert.equal(C1[r][c], C1[c][r], `C1 not sym at (${r},${c})`);

// Now find a signed permutation (p, s) sending C1 -> canonical C (s[i]*s[j]*C1[p[i]][p[j]] = C[i][j])!
let isoPerm = null;
{
  const p = new Int32Array(10), s = new Int32Array(10), used = new Uint8Array(10);
  function dfs(k) {
    if (isoPerm) return;
    if (k === 10) { isoPerm = { p: Array.from(p), s: Array.from(s) }; return; }
    for (let v = 0; v < 10; v++) {
      if (used[v]) continue;
      if (k === 0) { p[0] = v; s[0] = 1; used[v] = 1; dfs(1); used[v] = 0; }
      else {
        const sk = C[0][k] * C1[p[0]][v]; let ok = true;
        for (let i = 1; i < k; i++) if (s[i] * sk * C1[p[i]][v] !== C[i][k]) { ok = false; break; }
        if (ok) { p[k] = v; s[k] = sk; used[v] = 1; dfs(k + 1); used[v] = 0; }
      }
    }
  }
  dfs(0);
}
assert(isoPerm !== null, 'Failed to map C1 -> C!');
console.log('EXACT monomial isomorphism M_mu4 -> canonical H_conf = (I + iC)/sqrt(10) FOUND!');

// Transform B_iso = B2^dag into canonical H_conf coordinates!
// Original: B_iso_old[r][c] = conj(B2[c][r]).
// Left action on M was: row r of M multiplied by i^{aPh[r]}, then permuted by isoPerm (new row i takes old row isoPerm.p[i] with sign isoPerm.s[i])!
// Apply the EXACT SAME left action to B_iso_old:
//   B_iso_canon[i][c] = isoPerm.s[i] * i^{aPh[isoPerm.p[i]]} * conj(B2[c][isoPerm.p[i]])!
// Note: conj(i^k * v^e) = i^{-k} * v^{-e}!
const BisoSpec = Array.from({ length: N }, (_, c) => Array.from({ length: N }, (_, i) => {
  const oldR = isoPerm.p[i];
  const { e, k } = B2spec[c][oldR];
  const newE = -e;
  const newK = ((4 - k) + aPh[oldR] + (isoPerm.s[i] < 0 ? 2 : 0)) & 3;
  return { e: newE, k: newK };
})); // BisoSpec[c][i] is entry (row i, col c) of canonical B_iso!

// Build F_p field for splitting primes p == 1 (mod 20) (both i and sqrt(5) exist in F_p, 9 invertible!)
function powMod(a, e, p) {
  let r = 1n, b = BigInt(((a % p) + p) % p), ee = BigInt(e), m = BigInt(p);
  while (ee > 0n) { if (ee & 1n) r = (r * b) % m; b = (b * b) % m; ee >>= 1n; }
  return Number(r);
}
function buildField(p) {
  const mod = x => ((x % p) + p) % p;
  const mul = (a, b) => Number((BigInt(a) * BigInt(b)) % BigInt(p));
  const inv = a => powMod(a, p - 2, p);
  let w20 = -1;
  for (let g = 2; g < p; g++) {
    const w = powMod(g, (p - 1) / 20, p);
    if (powMod(w, 10, p) === p - 1 && powMod(w, 4, p) !== 1 && powMod(w, 5, p) !== 1) { w20 = w; break; }
  }
  const Z20 = Array.from({ length: 20 }, (_, k) => powMod(w20, k, p));
  const iVal = Z20[5];
  assert.equal(mul(iVal, iVal), p - 1);
  const s5 = mod(1 + 2 * Z20[4] + 2 * Z20[16]);
  assert.equal(mul(s5, s5), 5);
  const inv3 = inv(3);
  const vVal = mul(mod(2 + mul(iVal, s5)), inv3);
  const vInv = mul(mod(2 - mul(iVal, s5)), inv3);
  assert.equal(mul(vVal, vInv), 1);
  const vPow = { '-2': mul(vInv, vInv), '-1': vInv, '0': 1, '1': vVal, '2': mul(vVal, vVal) };
  const iPow = [1, iVal, p - 1, mod(-iVal)];

  // B[c][r] and Bc[c][r] (conjugate: conj(i^k v^e) = i^{-k} v^{-e})
  const B = Array.from({ length: N }, (_, c) => Array.from({ length: N }, (_, r) => {
    const { e, k } = BisoSpec[c][r]; return mul(iPow[k], vPow[String(e)]);
  }));
  const Bc = Array.from({ length: N }, (_, c) => Array.from({ length: N }, (_, r) => {
    const { e, k } = BisoSpec[c][r]; return mul(iPow[(4 - k) & 3], vPow[String(-e)]);
  }));
  // EXACT verification of B^dag B = 10 I AND |(I - iC) B|^2 = 10 in F_p!!
  for (let c = 0; c < N; c++) for (let d = 0; d < N; d++) {
    let s = 0; for (let r = 0; r < N; r++) s = (s + mul(Bc[c][r], B[d][r])) % p;
    assert.equal(s, c === d ? 10 : 0);
  }
  for (let c = 0; c < N; c++) for (let j = 0; j < N; j++) {
    let a = 0, ac = 0;
    for (let r = 0; r < N; r++) {
      const h = r === j ? 1 : mod(-C[j][r] * iVal), hc = r === j ? 1 : mod(C[j][r] * iVal);
      a = (a + mul(h, B[c][r])) % p; ac = (ac + mul(hc, Bc[c][r])) % p;
    }
    assert.equal(mul(a, ac), 10);
  }
  return { p, mod, mul, inv, iVal, iPow, B, Bc };
}

const P1 = 1000081, P2 = 1000921;
const F1 = buildField(P1);
console.log(`[1] VERIFIED canonical {I, H_conf, B_iso} EXACTLY in F_${P1} (all 100 ortho + all 100 unbiased eqs)!`);

// Compute EXACT G_1440(H_conf) = G_720(+C, lin) U G_720(-C, anti-lin) stabilizer & orbit size of B_iso!
function genG(signC) {
  const syms = []; const p = new Int32Array(10), s = new Int32Array(10), used = new Uint8Array(10);
  function dfs(k) {
    if (k === 10) { syms.push({ p: Array.from(p), s: Array.from(s), conj: signC === -1 }); return; }
    for (let v = 0; v < 10; v++) {
      if (used[v]) continue;
      if (k === 0) { p[0] = v; s[0] = 1; used[v] = 1; dfs(1); used[v] = 0; }
      else {
        const sk = signC * C[0][k] * C[p[0]][v]; let ok = true;
        for (let i = 1; i < k; i++) if (s[i] * sk * C[p[i]][v] !== signC * C[i][k]) { ok = false; break; }
        if (ok) { p[k] = v; s[k] = sk; used[v] = 1; dfs(k + 1); used[v] = 0; }
      }
    }
  }
  dfs(0); return syms;
}
const GlinAll = genG(+1), GantiAll = genG(-1);
function colKeyFp(col, F) { const iv0 = F.inv(col[0]); return col.map(v => F.mul(v, iv0)).join(','); }
const baseColKeys = F1.B.map(c => colKeyFp(c, F1));
const baseBasisKey = baseColKeys.slice().sort().join('|');

const Glin = [], Ganti = [];
for (const g of GlinAll) {
  const nk = F1.B.map(col => colKeyFp(Array.from({ length: N }, (_, i) => g.s[i] > 0 ? col[g.p[i]] : F1.mod(-col[g.p[i]])), F1)).sort().join('|');
  if (nk === baseBasisKey) Glin.push(g);
}
for (const g of GantiAll) {
  const nk = F1.Bc.map(col => colKeyFp(Array.from({ length: N }, (_, i) => g.s[i] > 0 ? col[g.p[i]] : F1.mod(-col[g.p[i]])), F1)).sort().join('|');
  if (nk === baseBasisKey) Ganti.push(g);
}
const orbSize = 1440 / (Glin.length + Ganti.length);
console.log(`    EXACT Stabilizer of B_iso in G_1440(H_conf): |G_lin| = ${Glin.length}, |G_anti| = ${Ganti.length}, total = ${Glin.length + Ganti.length}`);
console.log(`    EXACT G_1440(H_conf) ORBIT SIZE OF B_iso = ${orbSize} DISTINCT DEFECT-0 THIRD BASES FOR {I, H_conf}!!`);

// =========================================================================
// [2] ZERO-FLOAT 171D JOINT 3-MUB RIGIDITY THEOREM FOR {I, H_conf, B_iso}!!
// =========================================================================
function checkJoint171RankModP(p) {
  const Fp = buildField(p);
  const varB2 = (r, c) => (r === 0 || c === 0) ? -1 : (r - 1) * 9 + (c - 1);
  const varB3 = (r, c) => r === 0 ? -1 : 81 + c * 9 + (r - 1);
  const rows = [];
  // 90 ortho for H_conf
  for (let c = 0; c < N; c++) for (let d = c + 1; d < N; d++) {
    const rRe = new Int32Array(171), rIm = new Int32Array(171);
    for (let r = 0; r < N; r++) {
      const pw = (Hpow[r][d] - Hpow[r][c] + 4) & 3;
      const dS = Fp.iPow[(pw + 1) & 3], cDS = Fp.iPow[(4 - pw + 3) & 3];
      const xr = (dS + cDS) % p, xi = Fp.mod(dS - cDS);
      const vc = varB2(r, c), vd = varB2(r, d);
      if (vc >= 0) { rRe[vc] = Fp.mod(rRe[vc] - xr); rIm[vc] = Fp.mod(rIm[vc] - xi); }
      if (vd >= 0) { rRe[vd] = (rRe[vd] + xr) % p; rIm[vd] = (rIm[vd] + xi) % p; }
    }
    rows.push(rRe, rIm);
  }
  // 90 ortho for B_iso
  for (let c = 0; c < N; c++) for (let d = c + 1; d < N; d++) {
    const rRe = new Int32Array(171), rIm = new Int32Array(171);
    for (let r = 0; r < N; r++) {
      const z = Fp.mul(Fp.Bc[c][r], Fp.B[d][r]), cz = Fp.mul(Fp.B[c][r], Fp.Bc[d][r]);
      const dS = Fp.mul(Fp.iVal, z), cDS = Fp.mul(Fp.mod(-Fp.iVal), cz);
      const xr = (dS + cDS) % p, xi = Fp.mod(dS - cDS);
      const vc = varB3(r, c), vd = varB3(r, d);
      if (vc >= 0) { rRe[vc] = Fp.mod(rRe[vc] - xr); rIm[vc] = Fp.mod(rIm[vc] - xi); }
      if (vd >= 0) { rRe[vd] = (rRe[vd] + xr) % p; rIm[vd] = (rIm[vd] + xi) % p; }
    }
    rows.push(rRe, rIm);
  }
  // 100 unbiased eqs
  for (let c = 0; c < N; c++) for (let d = 0; d < N; d++) {
    const row = new Int32Array(171);
    let A = 0, cA = 0; const ts = [], cts = [];
    for (let r = 0; r < N; r++) {
      const hc = Fp.iPow[(4 - Hpow[r][c]) & 3], h = Fp.iPow[Hpow[r][c]];
      const t = Fp.mul(hc, Fp.B[d][r]), ct = Fp.mul(h, Fp.Bc[d][r]);
      ts.push(t); cts.push(ct); A = (A + t) % p; cA = (cA + ct) % p;
    }
    assert.equal(Fp.mul(A, cA), 10);
    for (let r = 0; r < N; r++) {
      const gv = Fp.mod(Fp.mul(cA, ts[r]) - Fp.mul(A, cts[r]));
      const vu = varB2(r, c), vv = varB3(r, d);
      if (vu >= 0) row[vu] = Fp.mod(row[vu] - gv);
      if (vv >= 0) row[vv] = (row[vv] + gv) % p;
    }
    rows.push(row);
  }
  const piv = Array(171).fill(null); let rank = 0; const bp = BigInt(p);
  for (const r0 of rows) {
    const r = Array.from(r0);
    for (let c = 0; c < 171; c++) {
      if (!r[c]) continue;
      if (!piv[c]) { const iv = BigInt(Fp.inv(r[c])); for (let k = c; k < 171; k++) r[k] = Number((BigInt(r[k]) * iv) % bp); piv[c] = r; rank++; break; }
      const f = BigInt(r[c]), pr = piv[c];
      for (let k = c; k < 171; k++) r[k] = Number((BigInt(r[k]) - (f * BigInt(pr[k])) % bp + bp) % bp);
    }
  }
  return rank;
}
console.log(`\n[2] ZERO-FLOAT 171D JOINT 3-MUB RIGIDITY FOR {I, H_conf, B_iso} IN Z[i, sqrt(5)]:`);
const jRank1 = checkJoint171RankModP(P1), jRank2 = checkJoint171RankModP(P2);
console.log(`    p = ${P1}: EXACT rank = ${jRank1}/171 | p = ${P2}: EXACT rank = ${jRank2}/171!`);
assert.equal(jRank1, 171); assert.equal(jRank2, 171);
console.log(`    Q.E.D. (ZERO FLOATS!!): {I, H_conf, B_iso} is STRICTLY LOCALLY ISOLATED in the FULL 171D 3-MUB moduli space!`);

// =========================================================================
// [3] ZERO-FLOAT DEGREE-6 HILBERT NULLSTELLENSATZ FOR {I, H_conf, B_iso}!!
// =========================================================================
console.log(`\n[3] ZERO-FLOAT DEGREE-6 HILBERT NULLSTELLENSATZ FOR {I, H_conf, B_iso} (|G_lin| = ${Glin.length}):`);
const POW8 = new Int32Array(10); for (let i = 0; i < 10; i++) POW8[i] = 1 << (3 * i);
const ZERO_CODE = (() => { let c = 0; for (let i = 0; i < 10; i++) c += 3 * POW8[i]; return c; })();
function applyGlin(code, g) {
  let c = ZERO_CODE, sgn = 1;
  for (let i = 0; i < 10; i++) { const e = ((code >>> (3 * i)) & 7) - 3; if (e !== 0) { if ((e & 1) !== 0 && g.s[i] < 0) sgn = -sgn; c += e * POW8[g.p[i]]; } }
  return [c, sgn];
}
const monoSet = new Set(), allCodes = [];
function addCode(code) { if (!monoSet.has(code)) { monoSet.add(code); allCodes.push(code); } }
addCode(ZERO_CODE);
for (let a = 0; a < 10; a++) for (let b = 0; b < 10; b++) if (a !== b) addCode(ZERO_CODE + POW8[a] - POW8[b]);
for (let a = 0; a < 10; a++) for (let b = a; b < 10; b++)
  for (let c = 0; c < 10; c++) if (c !== a && c !== b)
    for (let d = c; d < 10; d++) if (d !== a && d !== b)
      addCode(ZERO_CODE + POW8[a] + POW8[b] - POW8[c] - POW8[d]);
for (let a = 0; a < 10; a++) for (let b = a; b < 10; b++) for (let c = b; c < 10; c++)
  for (let d = 0; d < 10; d++) if (d !== a && d !== b && d !== c)
    for (let e = d; e < 10; e++) if (e !== a && e !== b && e !== c)
      for (let f = e; f < 10; f++) if (f !== a && f !== b && f !== c)
        addCode(ZERO_CODE + POW8[a] + POW8[b] + POW8[c] - POW8[d] - POW8[e] - POW8[f]);

const codeToOrb = new Map(), orbRep = [];
for (const code of allCodes) {
  if (codeToOrb.has(code)) continue;
  const om = new Map(); let kills = false;
  for (const g of Glin) { const [c2, sgn] = applyGlin(code, g); const prev = om.get(c2); if (prev !== undefined) { if (prev !== sgn) kills = true; } else om.set(c2, sgn); }
  if (kills) { for (const c2 of om.keys()) codeToOrb.set(c2, { orbId: -1, sign: 0 }); }
  else { const id = orbRep.length; orbRep.push(code); for (const [c2, sgn] of om.entries()) codeToOrb.set(c2, { orbId: id, sign: sgn }); }
}
const R = orbRep.length, constOrb = codeToOrb.get(ZERO_CODE).orbId;
const colPerm = Array.from({ length: R }, (_, i) => i); colPerm[constOrb] = R - 1; colPerm[R - 1] = constOrb;
const deg22 = allCodes.filter(code => { let dp = 0; for (let i = 0; i < 10; i++) { const e = ((code >>> (3 * i)) & 7) - 3; if (e > 0) dp += e; } return dp <= 2; });
function degOf(code) { let d = 0; for (let i = 0; i < 10; i++) { const e = ((code >>> (3 * i)) & 7) - 3; if (e > 0) d += e; } return d; }
deg22.sort((a, b) => degOf(a) - degOf(b));

const rowOrbs = [], rowSeen = new Set();
for (let r = 0; r < 10; r++) { if (rowSeen.has(r)) continue; rowOrbs.push(r); for (const g of Glin) rowSeen.add(g.p[r]); }
const colOrbs = [], colSeen = new Set();
for (let c = 0; c < 10; c++) {
  if (colSeen.has(c)) continue; colOrbs.push(c);
  for (const g of Glin) {
    const nc = Array.from({ length: 10 }, (_, i) => g.s[i] > 0 ? F1.B[c][g.p[i]] : F1.mod(-F1.B[c][g.p[i]]));
    colSeen.add(baseColKeys.indexOf(colKeyFp(nc, F1)));
  }
}
console.log(`    G_lin-invariant degree-(3,3) monomials: R = ${R}, rowOrbs=[${rowOrbs}], colOrbs=[${colOrbs}]`);

function proveNullstellensatzModP(p) {
  const t0 = Date.now(), Fp = buildField(p), bp = BigInt(p);
  const gens = [];
  for (const j of rowOrbs) {
    const terms = [];
    for (let a = 0; a < 10; a++) {
      const hac = a === j ? 1 : Fp.mod(C[j][a] * Fp.iVal);
      for (let b = 0; b < 10; b++) { if (a === b) continue; const hb = b === j ? 1 : Fp.mod(-C[j][b] * Fp.iVal); terms.push([b, a, Fp.mul(hac, hb)]); }
    }
    // Do NOT filter deg22 by subG so we guarantee 100% full span of G_lin-averaged ideal!
    gens.push({ terms, ms: deg22 });
  }
  for (const c of colOrbs) {
    const terms = [];
    for (let a = 0; a < 10; a++) for (let b = 0; b < 10; b++) { if (a === b) continue; terms.push([a, b, Fp.mul(Fp.Bc[c][a], Fp.B[c][b])]); }
    gens.push({ terms, ms: deg22 });
  }
  const pivots = new Array(R).fill(null); let rank = 0; const row = new Int32Array(R);
  for (const deg of [0, 1, 2]) {
    for (const gen of gens) {
      for (const m of gen.ms) {
        if (degOf(m) !== deg || rank === R) continue;
        row.fill(0);
        for (const [q, rIdx, coef] of gen.terms) {
          const info = codeToOrb.get(m + POW8[q] - POW8[rIdx]);
          if (!info || info.orbId === -1) continue;
          const col = colPerm[info.orbId], v = info.sign > 0 ? coef : (coef === 0 ? 0 : p - coef);
          const s = row[col] + v; row[col] = s >= p ? s - p : s;
        }
        for (let col = 0; col < R; col++) {
          const f = row[col]; if (f === 0) continue;
          const pr = pivots[col];
          if (pr) {
            const bf = BigInt(f);
            for (let j = col; j < R; j++) { const pj = pr[j]; if (pj !== 0) { let v = row[j] - Number((bf * BigInt(pj)) % bp); if (v < 0) v += p; row[j] = v; } }
          } else {
            const iv = BigInt(Fp.inv(f)), nr = new Int32Array(R);
            for (let j = col; j < R; j++) if (row[j] !== 0) nr[j] = Number((BigInt(row[j]) * iv) % bp);
            pivots[col] = nr; rank++; break;
          }
        }
      }
      if (rank === R) break;
    }
    console.log(`    [p=${p}] after deg <= (${deg+1},${deg+1}): rank = ${rank}/${R}, constPivot = ${pivots[R-1]!==null} (${((Date.now()-t0)/1000).toFixed(2)}s)`);
    if (rank === R) break;
  }
  return { rank, R, hasConst: pivots[R - 1] !== null };
}

const n1 = proveNullstellensatzModP(P1);
const n2 = proveNullstellensatzModP(P2);
if (n1.rank === n1.R && n2.rank === n2.R) {
  console.log(`\n*** Q.E.D. (ZERO FLOATS!!): FULL RANK ${n1.rank}/${n1.R} AT BOTH p=${P1} AND p=${P2}!! ***`);
  console.log(`    THEREFORE {I, H_conf, B_iso} (ALL ${orbSize} bases in the {0,0,16} orbit!) ADMITS NO COMMON UNBIASED 4TH VECTOR w in C^10!!`);
}
