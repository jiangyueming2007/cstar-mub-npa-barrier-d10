// Round 13 Step 2: 100% ZERO-FLOAT DEGREE-6 HILBERT NULLSTELLENSATZ PROOF FOR {I, H_conf, B_new} (Orbit #0)!
// Uses ONLY the purely LINEAR stabilizer G_lin (order 20, NO conjugation of K-coefficients!) so every
// G_lin-symmetrized ideal element S(m, f) = sum_{g in G_lin} m(P_g y) f(P_g y) is strictly a K-linear
// combination of ideal generators in K[y^{+-1}], K = Q(zeta_20, sqrt(6))!
// If rank_{F_p}(M) == R (FULL ROW RANK over ALL R G_lin-invariant degree-(3,3) Laurent monomials),
// then an R x R submatrix M_0 has det(M_0) != 0 (mod p) => det(M_0) != 0 in K subset C => 1 in Ideal(C), Q.E.D.!

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
const BSPEC_STR = [
  ' 0  4  7 +2 +8  3  6 +6 +5 +4',
  ' 0  6  3 +8 +2  7  4 +4 +5 +6',
  ' 0  4  7 +7 +3  3  6 +1 +0 +9',
  ' 0  6  3 +3 +7  7  4 +9 +0 +1',
  ' 0  8  9 -4 -6  1  2 -2 -5 -8',
  ' 0  8  9 -9 -1  1  2 -7 -0 -3',
  ' 0  2  1 -6 -4  9  8 -8 -5 -2',
  ' 0  2  1 -1 -9  9  8 -3 -0 -7',
  ' 0  0  5  0  0  5  0  0  5  0',
  ' 0  0  5  5  5  5  0  5  0  5',
];
const Bspec = BSPEC_STR.map(line => line.trim().split(/\s+/).map(tok => {
  if (tok.startsWith('+')) return { e: 1, m: Number(tok.slice(1)) };
  if (tok.startsWith('-')) return { e: -1, m: Number(tok.slice(1)) };
  return { e: 0, m: Number(tok) };
}));

function powMod(a, e, p) {
  let r = 1n, b = BigInt(((a % p) + p) % p), ee = BigInt(e), m = BigInt(p);
  while (ee > 0n) { if (ee & 1n) r = (r * b) % m; b = (b * b) % m; ee >>= 1n; }
  return Number(r);
}
function sqrtModP(a, p) {
  a = ((a % p) + p) % p;
  let q = p - 1, s = 0; while ((q & 1) === 0) { q >>= 1; s++; }
  let z = 2; while (powMod(z, (p - 1) / 2, p) !== p - 1) z++;
  let m = s, c = powMod(z, q, p), t = powMod(a, q, p), r = powMod(a, (q + 1) / 2, p);
  while (t !== 1) {
    let i = 1, tt = Number((BigInt(t) * BigInt(t)) % BigInt(p));
    while (tt !== 1) { tt = Number((BigInt(tt) * BigInt(tt)) % BigInt(p)); i++; }
    const b = powMod(c, 1 << (m - i - 1), p);
    m = i; c = Number((BigInt(b) * BigInt(b)) % BigInt(p));
    t = Number((BigInt(t) * BigInt(c)) % BigInt(p));
    r = Number((BigInt(r) * BigInt(b)) % BigInt(p));
  }
  assert.equal(Number((BigInt(r) * BigInt(r)) % BigInt(p)), a);
  return r;
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
  const s6 = sqrtModP(6, p), invS6 = inv(s6);
  const uVal = mul(mod(1 + mul(iVal, s5)), invS6);
  const uInv = mul(mod(1 - mul(iVal, s5)), invS6);
  assert.equal(mul(uVal, uInv), 1);

  const B = Array.from({ length: 10 }, (_, c) => Array.from({ length: 10 }, (_, r) => {
    const { e, m } = Bspec[c][r], z = Z20[(2 * m) % 20];
    return e === 0 ? z : mul(e === 1 ? uVal : uInv, z);
  }));
  const Bc = Array.from({ length: 10 }, (_, c) => Array.from({ length: 10 }, (_, r) => {
    const { e, m } = Bspec[c][r], z = Z20[(20 - 2 * m) % 20];
    return e === 0 ? z : mul(e === 1 ? uInv : uVal, z);
  }));
  for (let c = 0; c < 10; c++) for (let d = 0; d < 10; d++) {
    let s = 0; for (let r = 0; r < 10; r++) s = (s + mul(Bc[c][r], B[d][r])) % p;
    assert.equal(s, c === d ? 10 : 0);
  }
  for (let c = 0; c < 10; c++) for (let j = 0; j < 10; j++) {
    let a = 0, ac = 0;
    for (let r = 0; r < 10; r++) {
      const h = r === j ? 1 : mod(-C[j][r] * iVal), hc = r === j ? 1 : mod(C[j][r] * iVal);
      a = (a + mul(h, B[c][r])) % p; ac = (ac + mul(hc, Bc[c][r])) % p;
    }
    assert.equal(mul(a, ac), 10);
  }
  return { p, mod, mul, inv, iVal, B, Bc };
}

const P1 = 1000081, P2 = 1000921;
const F1 = buildField(P1);

// 1. Build ONLY the purely LINEAR stabilizer G_lin (order 20, NO complex conjugation!)
function colKeyFp(col, F) { const iv0 = F.inv(col[0]); return col.map(v => F.mul(v, iv0)).join(','); }
const baseColKeys = F1.B.map(col => colKeyFp(col, F1));
const baseBasisKey = baseColKeys.slice().sort().join('|');
const Glin = [];
const pp = new Int32Array(10), ss = new Int32Array(10), used = new Uint8Array(10);
(function dfsLin(k) {
  if (k === 10) {
    const nKeys = F1.B.map(col => colKeyFp(Array.from({ length: 10 }, (_, i) => ss[i] > 0 ? col[pp[i]] : F1.mod(-col[pp[i]])), F1)).sort().join('|');
    if (nKeys === baseBasisKey) Glin.push({ p: Array.from(pp), s: Array.from(ss) });
    return;
  }
  for (let v = 0; v < 10; v++) {
    if (used[v]) continue;
    if (k === 0) { pp[0] = v; ss[0] = 1; used[v] = 1; dfsLin(1); used[v] = 0; }
    else {
      const sk = C[pp[0]][v]; let ok = true;
      for (let i = 1; i < k; i++) if (ss[i] * sk * C[pp[i]][v] !== C[i][k]) { ok = false; break; }
      if (ok) { pp[k] = v; ss[k] = sk; used[v] = 1; dfsLin(k + 1); used[v] = 0; }
    }
  }
})(0);
assert.equal(Glin.length, 20);
console.log(`Purely linear stabilizer |G_lin| = ${Glin.length} (exact order 20, ZERO conjugation!).`);

// Find generator orbits under G_lin:
// Note: (P_g y)_i = s[i] * y_{p[i]}, so for a vector v, <v, P_g y> = sum_i conj(v_i) s[i] y_{p[i]} = <P_g^T v, y>
// Wait: C_{p[i], p[j]} = s[i]*s[j]*C_{i,j} means (P_g C P_g^T)_{a,b} = C_{a,b}, so P_g preserves C and H_conf!
const rowOrbs = [], rowSeen = new Set();
for (let r = 0; r < 10; r++) {
  if (rowSeen.has(r)) continue;
  rowOrbs.push(r);
  for (const g of Glin) rowSeen.add(g.p[r]);
}
const colOrbs = [], colSeen = new Set();
for (let c = 0; c < 10; c++) {
  if (colSeen.has(c)) continue;
  colOrbs.push(c);
  for (const g of Glin) {
    const nc = Array.from({ length: 10 }, (_, i) => g.s[i] > 0 ? F1.B[c][g.p[i]] : F1.mod(-F1.B[c][g.p[i]]));
    colSeen.add(baseColKeys.indexOf(colKeyFp(nc, F1)));
  }
}
console.log(`G_lin generator orbits: row reps = [${rowOrbs}], col reps = [${colOrbs}]`);

// 2. Collapse all 24,691 degree <= (3,3) Laurent monomials under G_lin!
const POW8 = new Int32Array(10); for (let i = 0; i < 10; i++) POW8[i] = 1 << (3 * i);
const ZERO_CODE = (() => { let c = 0; for (let i = 0; i < 10; i++) c += 3 * POW8[i]; return c; })();
function applyGlin(code, g) {
  let c = ZERO_CODE, sgn = 1;
  for (let i = 0; i < 10; i++) {
    const e = ((code >>> (3 * i)) & 7) - 3;
    if (e !== 0) { if ((e & 1) !== 0 && g.s[i] < 0) sgn = -sgn; c += e * POW8[g.p[i]]; }
  }
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
  for (const g of Glin) {
    const [c2, sgn] = applyGlin(code, g);
    const prev = om.get(c2);
    if (prev !== undefined) { if (prev !== sgn) kills = true; } else om.set(c2, sgn);
  }
  if (kills) { for (const c2 of om.keys()) codeToOrb.set(c2, { orbId: -1, sign: 0 }); }
  else { const id = orbRep.length; orbRep.push(code); for (const [c2, sgn] of om.entries()) codeToOrb.set(c2, { orbId: id, sign: sgn }); }
}
const R = orbRep.length;
const constOrb = codeToOrb.get(ZERO_CODE).orbId;
console.log(`G_lin-invariant degree-(3,3) Laurent monomials: R = ${R} (constOrb = ${constOrb})`);

// Order columns so constOrb is LAST (index R-1) — then checking BOTH full rank R/R AND whether e_const alone is in span is trivial!
const colPerm = Array.from({ length: R }, (_, i) => i);
colPerm[constOrb] = R - 1; colPerm[R - 1] = constOrb;
const deg22 = [];
for (const code of allCodes) { let dp = 0; for (let i = 0; i < 10; i++) { const e = ((code >>> (3 * i)) & 7) - 3; if (e > 0) dp += e; } if (dp <= 2) deg22.push(code); }
function degOf(code) { let d = 0; for (let i = 0; i < 10; i++) { const e = ((code >>> (3 * i)) & 7) - 3; if (e > 0) d += e; } return d; }
deg22.sort((a, b) => degOf(a) - degOf(b));

function proveModP(p) {
  const t0 = Date.now();
  const Fp = buildField(p);
  const gens = [];
  for (const j of rowOrbs) {
    const terms = [];
    for (let a = 0; a < 10; a++) {
      const hac = a === j ? 1 : Fp.mod(C[j][a] * Fp.iVal);
      for (let b = 0; b < 10; b++) {
        if (a === b) continue;
        const hb = b === j ? 1 : Fp.mod(-C[j][b] * Fp.iVal);
        terms.push([b, a, Fp.mul(hac, hb)]);
      }
    }
    // Filter multipliers by stabilizer of row j in G_lin (optional speedup, exact)
    const subG = Glin.filter(g => g.p[j] === j);
    const seen = new Set(), ms = [];
    for (const m of deg22) { if (seen.has(m)) continue; ms.push(m); for (const g of subG) seen.add(applyGlin(m, g)[0]); }
    gens.push({ name: `h_${j}`, terms, ms });
  }
  for (const c of colOrbs) {
    const terms = [];
    for (let a = 0; a < 10; a++) for (let b = 0; b < 10; b++) {
      if (a === b) continue;
      terms.push([a, b, Fp.mul(Fp.Bc[c][a], Fp.B[c][b])]);
    }
    const cKey = colKeyFp(Fp.B[c], Fp);
    const subG = Glin.filter(g => {
      const nc = Array.from({ length: 10 }, (_, i) => g.s[i] > 0 ? Fp.B[c][g.p[i]] : Fp.mod(-Fp.B[c][g.p[i]]));
      return colKeyFp(nc, Fp) === cKey;
    });
    const seen = new Set(), ms = [];
    for (const m of deg22) { if (seen.has(m)) continue; ms.push(m); for (const g of subG) seen.add(applyGlin(m, g)[0]); }
    gens.push({ name: `b_${c}`, terms, ms });
  }

  const pivots = new Array(R).fill(null);
  let rank = 0;
  const bp = BigInt(p);
  const row = new Int32Array(R);

  for (const deg of [0, 1, 2]) {
    for (const gen of gens) {
      for (const m of gen.ms) {
        if (degOf(m) !== deg) continue;
        if (rank === R) break;
        row.fill(0);
        for (const [q, rIdx, coef] of gen.terms) {
          const info = codeToOrb.get(m + POW8[q] - POW8[rIdx]);
          if (!info || info.orbId === -1) continue;
          const col = colPerm[info.orbId];
          const v = info.sign > 0 ? coef : (coef === 0 ? 0 : p - coef);
          const s = row[col] + v;
          row[col] = s >= p ? s - p : s;
        }
        // Reduce row
        for (let col = 0; col < R; col++) {
          const f = row[col]; if (f === 0) continue;
          const pr = pivots[col];
          if (pr) {
            const bf = BigInt(f);
            for (let j = col; j < R; j++) {
              const pj = pr[j];
              if (pj !== 0) {
                let v = row[j] - Number((bf * BigInt(pj)) % bp);
                if (v < 0) v += p;
                row[j] = v;
              }
            }
          } else {
            const iv = BigInt(Fp.inv(f));
            const nr = new Int32Array(R);
            for (let j = col; j < R; j++) if (row[j] !== 0) nr[j] = Number((BigInt(row[j]) * iv) % bp);
            pivots[col] = nr;
            rank++;
            break;
          }
        }
      }
      if (rank === R) break;
    }
    console.log(`    [p=${p}] after degree <= (${deg+1},${deg+1}): rank = ${rank} / ${R}, constPivot = ${pivots[R - 1] !== null} (${((Date.now() - t0)/1000).toFixed(2)}s)`);
    if (rank === R) break;
  }
  return { p, rank, R, hasConstPivot: pivots[R - 1] !== null, ms: Date.now() - t0 };
}

console.log('\nRunning ZERO-FLOAT Degree-6 Hilbert Nullstellensatz proof for {I, H_conf, B_new} (Orbit #0):');
const res1 = proveModP(P1);
const res2 = proveModP(P2);
if (res1.rank === res1.R && res2.rank === res2.R) {
  console.log(`\n*** Q.E.D. (ZERO FLOATS!!): FULL RANK ${res1.rank}/${res1.R} AT BOTH p=${P1} AND p=${P2}!! ***`);
  console.log(`    An ${res1.R}x${res1.R} submatrix M_0 of degree-(3,3) ideal generators over K=Q(zeta_20,sqrt(6))`);
  console.log(`    has det(M_0) != 0 (mod p), hence M_0 is INVERTIBLE over C => 1 in <h_0..h_9, b_0..b_9>_{<=6} over C!`);
  console.log(`    THEREFORE {I, H_conf, B_new} (ALL 72 bases in Orbit #0) ADMITS NO COMMON UNBIASED 4TH VECTOR w in C^10!!`);
} else if (res1.hasConstPivot) {
  console.log(`\n*** e_const IS IN THE DEGREE-6 IDEAL SPAN (pivot at constOrb = R-1, rank = ${res1.rank}/${res1.R})!! ***`);
} else {
  console.log(`\nDegree-6 rank = ${res1.rank}/${res1.R}, constPivot = ${res1.hasConstPivot}.`);
}
