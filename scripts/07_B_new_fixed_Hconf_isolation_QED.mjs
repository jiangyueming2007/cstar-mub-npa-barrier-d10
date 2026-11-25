import _pkgPath from 'node:path';
import { fileURLToPath as _pkgFileURLToPath } from 'node:url';
const _pkgDir = _pkgPath.dirname(_pkgFileURLToPath(import.meta.url));
// Round 20 Prong 2: ZERO-FLOAT PROOF OF STRICT FIXED-{I, H_conf} ISOLATION OF B_new (90/90 FULL RANK)
// AND STRICT TRANSVERSALITY OF THE 4D MANIFOLD M_4!
//
// Clarifies a subtle point raised in peer review:
//   Although {I, H_conf, B_new} has Joint 171D 3-MUB defect = 4 (lying on the 4D continuous manifold M_4),
//   when the pair {I, H_conf} is HELD FIXED, is B_new part of a continuous family of third bases, or STRICTLY ISOLATED?
//   Here we compute the exact 90-variable fixed-pair Jacobian (190 eqs x 90 vars) and 171-variable joint Jacobian
//   (280 eqs x 171 vars) modulo TWO independent split primes p1, p2 (where 20 | (p-1) and 6, 10 are quadratic residues)
//   with ZERO floating-point arithmetic, proving:
//     1. rank_{F_p}(J_fixed(B_new)) == 90 / 90  =>  B_new is a STRICTLY ISOLATED (0-dimensional) third basis of {I, H_conf}!
//     2. rank_{F_p}(J_joint({I, H_conf, B_new})) == 167 / 171  =>  exact joint defect <= 4 (and = 4 via Round 12/19 smooth curves).
//     3. Every non-zero tangent vector in ker(J_joint) has strictly non-zero B_2-component (purely transverse to fixed H_conf!).

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

const Bdata = JSON.parse(fs.readFileSync(_pkgPath.join(_pkgDir, '..', 'data', 'B_new_exotic_third_basis.json'), 'utf8'));
const Bcols = Bdata.columnsAnglesOverPi;

// Exact symbolic extraction B_new[r,c] = u^e * zeta_10^m, u = (1 + i*sqrt(5))/sqrt(6)
const argU = Math.atan2(Math.sqrt(5), 1);
const spec = [];
for (let c = 0; c < 10; c++) {
  const col = [];
  for (let r = 0; r < 10; r++) {
    const th = Bcols[c][r] * Math.PI;
    let best = null, bestErr = 1e9;
    for (const e of [-1, 0, 1]) for (let m = 0; m < 10; m++) {
      const cand = e * argU + m * (Math.PI / 5);
      const err = Math.hypot(Math.cos(th) - Math.cos(cand), Math.sin(th) - Math.sin(cand));
      if (err < bestErr) { bestErr = err; best = { e, m }; }
    }
    assert(bestErr < 1e-12);
    col.push(best);
  }
  spec.push(col);
}

// Modular arithmetic helpers (BigInt)
const mod = (a, p) => ((a % p) + p) % p;
function powMod(a, e, p) {
  let r = 1n, b = mod(a, p);
  while (e > 0n) { if (e & 1n) r = (r * b) % p; b = (b * b) % p; e >>= 1n; }
  return r;
}
const invMod = (a, p) => powMod(a, p - 2n, p);

function isPrime(n) {
  if (n < 2n) return false;
  for (let d = 2n; d * d <= n; d++) if (n % d === 0n) return false;
  return true;
}
function sqrtMod(a, p) {
  a = mod(a, p);
  if (a === 0n) return 0n;
  if (powMod(a, (p - 1n) / 2n, p) !== 1n) return null;
  // Tonelli-Shanks
  let q = p - 1n, s = 0n;
  while ((q & 1n) === 0n) { q >>= 1n; s++; }
  let z = 2n;
  while (powMod(z, (p - 1n) / 2n, p) !== p - 1n) z++;
  let m = s, c = powMod(z, q, p), t = powMod(a, q, p), r = powMod(a, (q + 1n) / 2n, p);
  while (t !== 1n) {
    let i = 1n, tmp = (t * t) % p;
    while (tmp !== 1n) { tmp = (tmp * tmp) % p; i++; }
    const b = powMod(c, 1n << (m - i - 1n), p);
    m = i; c = (b * b) % p; t = (t * c) % p; r = (r * b) % p;
  }
  assert((r * r) % p === a);
  return r;
}

// Find two split primes > 1,000,000 where 20 | (p-1) AND sqrt(5), sqrt(6), sqrt(10) all exist in F_p!
// Note: if 20 | (p-1), then -1 and 5 are automatically squares mod p! We just need 6 and 10 (i.e. 2 and 3) to be squares mod p!
const splitPrimes = [];
for (let cand = 1000001n; splitPrimes.length < 2; cand += 20n) {
  if (!isPrime(cand)) continue;
  if (sqrtMod(6n, cand) !== null && sqrtMod(10n, cand) !== null && sqrtMod(5n, cand) !== null && sqrtMod(cand - 1n, cand) !== null) {
    splitPrimes.push(cand);
  }
}

// Exact rank over F_p via Gaussian elimination, returning { rank, nullspace }
function rrefRankAndNull(rows, nCols, p) {
  const M = rows.map(r => r.slice());
  const nRows = M.length;
  const pivCols = [];
  let r = 0;
  for (let c = 0; c < nCols && r < nRows; c++) {
    let piv = -1;
    for (let i = r; i < nRows; i++) if (M[i][c] !== 0n) { piv = i; break; }
    if (piv === -1) continue;
    [M[r], M[piv]] = [M[piv], M[r]];
    const inv = invMod(M[r][c], p);
    for (let j = c; j < nCols; j++) M[r][j] = (M[r][j] * inv) % p;
    for (let i = 0; i < nRows; i++) {
      if (i === r || M[i][c] === 0n) continue;
      const f = M[i][c];
      for (let j = c; j < nCols; j++) M[i][j] = mod(M[i][j] - f * M[r][j], p);
    }
    pivCols.push(c);
    r++;
  }
  const pivSet = new Set(pivCols);
  const freeCols = [];
  for (let c = 0; c < nCols; c++) if (!pivSet.has(c)) freeCols.push(c);
  const nullspace = freeCols.map(fc => {
    const v = new Array(nCols).fill(0n);
    v[fc] = 1n;
    for (let i = 0; i < pivCols.length; i++) v[pivCols[i]] = mod(-M[i][fc], p);
    return v;
  });
  return { rank: r, nullspace };
}

function auditAtPrime(p) {
  // Find primitive 10th root of unity z10 in F_p, plus i_fp (i^2=-1), s5 (s5^2=5), s6 (s6^2=6), s10 (s10^2=10)
  const iFP = sqrtMod(p - 1n, p);
  const s6 = sqrtMod(6n, p);
  const s10 = sqrtMod(10n, p);
  // Find primitive 5th root of unity z5 in F_p -> z10 = -z5^3, and s5 = 1 + 2*z5 + 2*z5^4 (matching exact algebraic embedding!)
  let z5 = null;
  for (let g = 2n; g < p; g++) {
    const cand = powMod(g, (p - 1n) / 5n, p);
    if (cand !== 1n) { z5 = cand; break; }
  }
  assert(z5 !== null && powMod(z5, 5n, p) === 1n);
  const z5_4 = powMod(z5, 4n, p);
  const s5 = mod(1n + 2n * z5 + 2n * z5_4, p);
  assert((s5 * s5) % p === 5n);
  const z10_1 = mod(-powMod(z5, 3n, p), p);
  const Z10 = [1n];
  for (let m = 1; m < 10; m++) Z10.push((Z10[m - 1] * z10_1) % p);
  assert((Z10[9] * z10_1) % p === 1n);

  // Split real/imag parts of B_new[r, c] = cos(th) + i*sin(th) over F_p!
  // Remember: for any unimodular complex number z = x + i*y with conj(z) = z^{-1} = x - i*y:
  //   x = (z + z^{-1}) / 2,  y = (z - z^{-1}) / (2*i)!
  const inv2 = invMod(2n, p);
  const inv2i = invMod(mod(2n * iFP, p), p);
  const uFP = mod((1n + iFP * s5) * invMod(s6, p), p);
  const uInvFP = mod((1n - iFP * s5) * invMod(s6, p), p);
  assert((uFP * uInvFP) % p === 1n);

  const br = Array.from({ length: 10 }, () => new Array(10).fill(0n));
  const bi = Array.from({ length: 10 }, () => new Array(10).fill(0n));
  for (let c = 0; c < 10; c++) for (let r = 0; r < 10; r++) {
    const { e, m } = spec[c][r];
    const z = mod((e === 0 ? 1n : (e === 1 ? uFP : uInvFP)) * Z10[m], p);
    const zBar = mod((e === 0 ? 1n : (e === 1 ? uInvFP : uFP)) * Z10[(10 - m) % 10], p);
    assert((z * zBar) % p === 1n);
    br[c][r] = mod((z + zBar) * inv2, p);
    bi[c][r] = mod((z - zBar) * inv2i, p);
    assert(mod(br[c][r] * br[c][r] + bi[c][r] * bi[c][r], p) === 1n);
  }

  // H_conf[r, c] = (delta_{rc} + i*C_{rc}) / sqrt(10) -> multiply by sqrt(10) to get unimodular phases?
  // Wait: H_conf is NOT dephased (diag is 1/sqrt(10), off-diag is +-i/sqrt(10)), wait: |H_conf[r,c]| = 1/sqrt(10) for ALL r, c!
  // So sqrt(10) * H_conf[r, c] IS unimodular: hr[c][r] = (r===c ? 1 : 0), hi[c][r] = C[r][c] (already in {-1, 0, 1}!)
  // Wait: dephase H_conf so row 0 and col 0 are all 1!
  // Row 0 of sqrt(10)*H_conf is [1, i, i, ..., i]; Col 0 is [1, i, i, ..., i]^T!
  // Multiply row r>=1 by (-i) and col c>=1 by (-i) -> H_dep[r, c]:
  //   r=0 or c=0: 1
  //   r>=1, c>=1: (-i)*(-i)*(delta_{rc} + i*C_{rc}) = -(delta_{rc} + i*C_{rc})!
  // Apply the SAME left diagonal (-i for r>=1) to B_new so {I, H_dep, B_dep} is an exact dephased 3-MUB!
  const hr = Array.from({ length: 10 }, () => new Array(10).fill(0n));
  const hi = Array.from({ length: 10 }, () => new Array(10).fill(0n));
  const bdr = Array.from({ length: 10 }, () => new Array(10).fill(0n));
  const bdi = Array.from({ length: 10 }, () => new Array(10).fill(0n));
  for (let c = 0; c < 10; c++) for (let r = 0; r < 10; r++) {
    if (r === 0 || c === 0) { hr[c][r] = 1n; hi[c][r] = 0n; }
    else { hr[c][r] = (r === c) ? mod(-1n, p) : 0n; hi[c][r] = mod(-BigInt(C[r][c]), p); }
    // Left multiply B_new row r>=1 by (-i): (x + iy)*(-i) = y - ix
    if (r === 0) { bdr[c][r] = br[c][r]; bdi[c][r] = bi[c][r]; }
    else { bdr[c][r] = bi[c][r]; bdi[c][r] = mod(-br[c][r], p); }
  }
  // Verify exact orthogonality of H_dep, B_dep, and mutual unbiasedness |<H_dep[:,c], B_dep[:,d]>|^2 = 10 mod p!
  for (let c = 0; c < 10; c++) for (let d = c + 1; d < 10; d++) {
    let sr = 0n, si = 0n;
    for (let r = 0; r < 10; r++) {
      sr = mod(sr + bdr[c][r] * bdr[d][r] + bdi[c][r] * bdi[d][r], p);
      si = mod(si + bdr[c][r] * bdi[d][r] - bdi[c][r] * bdr[d][r], p);
    }
    assert(sr === 0n && si === 0n, `B_dep ortho failed`);
  }
  for (let c = 0; c < 10; c++) for (let d = 0; d < 10; d++) {
    let ar = 0n, ai = 0n;
    for (let r = 0; r < 10; r++) {
      ar = mod(ar + hr[c][r] * bdr[d][r] + hi[c][r] * bdi[d][r], p);
      ai = mod(ai + hr[c][r] * bdi[d][r] - hi[c][r] * bdr[d][r], p);
    }
    assert(mod(ar * ar + ai * ai, p) === 10n, `Unbiasedness failed`);
  }

  // Build EXACT 171-variable Joint 3-MUB Jacobian (280 rows x 171 cols):
  //   vars 0..80  : B2 = H_dep phases (r=1..9, c=1..9)
  //   vars 81..170: B3 = B_dep phases (c=0..9, r=1..9)
  // And the FIXED-{I, H_conf} 90-variable Jacobian is simply the 190 rows (B3 ortho + H-B3 unb) restricted to cols 81..170!
  const varB2 = (r, c) => (r === 0 || c === 0) ? -1 : (r - 1) * 9 + (c - 1);
  const varB3 = (r, c) => (r === 0) ? -1 : 81 + c * 9 + (r - 1);
  const jointRows = [];
  const fixedRows = [];

  const addOrthoRows = (ur, ui, vFn, targetArr, colOffset = 0, nCols = 171) => {
    for (let c = 0; c < 10; c++) for (let d = c + 1; d < 10; d++) {
      const rowRe = new Array(nCols).fill(0n), rowIm = new Array(nCols).fill(0n);
      for (let r = 0; r < 10; r++) {
        const pr = mod(ur[c][r] * ur[d][r] + ui[c][r] * ui[d][r], p);
        const pi = mod(ur[c][r] * ui[d][r] - ui[c][r] * ur[d][r], p);
        const vc = vFn(r, c), vd = vFn(r, d);
        if (vc >= 0) { rowRe[vc - colOffset] = mod(rowRe[vc - colOffset] + pi, p); rowIm[vc - colOffset] = mod(rowIm[vc - colOffset] - pr, p); }
        if (vd >= 0) { rowRe[vd - colOffset] = mod(rowRe[vd - colOffset] - pi, p); rowIm[vd - colOffset] = mod(rowIm[vd - colOffset] + pr, p); }
      }
      targetArr.push(rowRe, rowIm);
    }
  };
  const addUnbRows = (ur, ui, vU, vr, vi, vV, targetArr, colOffset = 0, nCols = 171) => {
    for (let c = 0; c < 10; c++) for (let d = 0; d < 10; d++) {
      let ar = 0n, ai = 0n;
      const tR = new Array(10), tI = new Array(10);
      for (let r = 0; r < 10; r++) {
        tR[r] = mod(ur[c][r] * vr[d][r] + ui[c][r] * vi[d][r], p);
        tI[r] = mod(ur[c][r] * vi[d][r] - ui[c][r] * vr[d][r], p);
        ar = mod(ar + tR[r], p); ai = mod(ai + tI[r], p);
      }
      const row = new Array(nCols).fill(0n);
      for (let r = 0; r < 10; r++) {
        const gv = mod(2n * (ai * tR[r] - ar * tI[r]), p);
        const vu = vU ? vU(r, c) : -1;
        const vv = vV(r, d);
        if (vu >= 0) row[vu - colOffset] = mod(row[vu - colOffset] - gv, p);
        if (vv >= 0) row[vv - colOffset] = mod(row[vv - colOffset] + gv, p);
      }
      targetArr.push(row);
    }
  };

  addOrthoRows(hr, hi, varB2, jointRows, 0, 171);
  addOrthoRows(bdr, bdi, varB3, jointRows, 0, 171);
  addUnbRows(hr, hi, varB2, bdr, bdi, varB3, jointRows, 0, 171);

  addOrthoRows(bdr, bdi, varB3, fixedRows, 81, 90);
  addUnbRows(hr, hi, null, bdr, bdi, varB3, fixedRows, 81, 90);

  const resFixed = rrefRankAndNull(fixedRows, 90, p);
  const resJoint = rrefRankAndNull(jointRows, 171, p);
  assert(resFixed.rank === 90, `Expected fixed 90D rank = 90, got ${resFixed.rank}`);
  assert(resJoint.rank === 167, `Expected joint 171D rank = 167, got ${resJoint.rank}`);

  // Verify strict transversality: EVERY vector in ker(J_joint) has nonzero projection onto B2 (cols 0..80),
  // in fact the 4 vectors projected to B2 are LINEARLY INDEPENDENT (rank 4/4)!
  const b2Projs = resJoint.nullspace.map(v => v.slice(0, 81));
  const b2ProjRank = rrefRankAndNull(b2Projs, 81, p).rank;
  assert(b2ProjRank === 4, `Expected B2-projection rank of M_4 tangent space = 4, got ${b2ProjRank}`);

  return {
    prime: p.toString(),
    fixed_Hconf_90D_rank: `${resFixed.rank}/90`,
    fixed_Hconf_Bnew_strictly_isolated: true,
    joint_171D_rank: `${resJoint.rank}/171`,
    joint_171D_defect: 171 - resJoint.rank,
    M4_tangent_B2_projection_rank: `${b2ProjRank}/4`,
    M4_strictly_transverse_to_fixed_Hconf: true,
  };
}

const t0 = Date.now();
console.log('=== ROUND 20: ZERO-FLOAT PROOF OF STRICT FIXED-{I, H_conf} ISOLATION OF B_new (90/90) & M_4 TRANSVERSALITY ===');
const out = splitPrimes.map(p => {
  const r = auditAtPrime(p);
  console.log(`  p = ${r.prime}:`, r);
  return r;
});
const elapsed = +((Date.now() - t0) / 1000).toFixed(2);
console.log(`Completed zero-float fixed-pair isolation & transversality proof in ${elapsed}s!`);
fs.writeFileSync(_pkgPath.join(_pkgDir, '..', 'certificates', '07_B_new_fixed_Hconf_isolation_QED.json'), JSON.stringify({ elapsedSec: elapsed, results: out }, null, 2));
