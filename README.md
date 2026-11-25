# Zero-Float Verification Package: $C^*$-Algebra MUBs ($N_{C^*}(d) = \infty$), Gribling–Polak Alignment, and Exact $d=10$ Impossibility Theorems

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24530887.svg)](https://doi.org/10.5281/zenodo.24530887)

**Author:** Yinghao Rui (芮颖浩)  
**Date:** November 25, 2026  
**Archival DOI (v1.0.0):** [10.5281/zenodo.24530887](https://doi.org/10.5281/zenodo.24530887) (All versions: [10.5281/zenodo.24530886](https://doi.org/10.5281/zenodo.24530886))  
**Repository:** https://github.com/jiangyueming2007/cstar-mub-npa-barrier-d10

This self-contained, **zero-dependency** verification package reproduces every zero-floating-point algebraic theorem in our two companion preprints:

1. **Paper II (Flagship arXiv-Ready LaTeX Preprint)**:  
   `Higher-Rank_Cstar-Algebra_MUBs_N_Cstar_Inf_Gribling-Polak_Alignment_and_Degree-2_Rank-1_Separation.tex`  
   *Higher-Rank $C^*$-Algebra Mutually Unbiased Bases ($N_{C^*}(d) = \infty$), Axiom-by-Axiom Gribling–Polak Alignment, and Universal Degree-2 Rank-1 Separation*
2. **Paper I (Companion English Preprint & Chinese Master Monograph)**:  
   - `Zero-Float_Impossibility_Theorems_for_d10_MUBs_and_Higher-Rank_Cstar_Hierarchy.md`
   - `十维互无偏基零浮点不可行定理与有限维C星代数MUB无穷阶阻断总论_芮颖浩.md`

---

## Requirements

- **Node.js v20+** (tested on Node.js v24.15.0)
- **Zero external packages** (`npm install` is NOT required; pure built-in `BigInt` + `node:assert` + `node:crypto`)
- Runs on any standard CPU (Windows, Linux, macOS) in **~2 minutes total**

---

## Quick Start (One-Click Verification)

```bash
# Run the complete 8-proof zero-float verification suite (~2 min)
node verify_all.mjs

# Run only Paper II (C*-MUB / NPA Barrier) proofs (~1.5 min)
node verify_all.mjs --paper2

# Run fast core proofs excluding d=22 (M_704) Clifford check (~35 sec)
node verify_all.mjs --fast
```

---

## Theorem-to-Script Mapping

| ID | Script | Theorem / Audit Verified (100% Zero-Float `BigInt`) | Approx. Time |
|---:|---|---|---:|
| `01` | `scripts/01_quaternionic_4_5_mub_QED.mjs` | **Tier 1 (Rank 2):** Closed-form Quaternionic 4-MUB ($d=6$) & 5-MUB ($d=10, 14$) in $\mathbb{Z}[i,\omega_p][j]$ + non-zero corner commutator | ~1s |
| `02` | `scripts/02_complex_M12_M20_projectors_QED.mjs` | Explicit complex rank-2 projectors in $M_{12}(\mathbb{Z}[i,\omega_3])$ (432 triple checks) & $M_{20}(\mathbb{Z}[i,\omega_5])$ (2,000 triple checks) | ~13s |
| `03` | `scripts/03_gribling_polak_def31_alignment_QED.mjs` | **Gribling–Polak (*Quantum* 2024, Def. 3.1) Axiom-by-Axiom Audit:** 144/500 intra-basis ortho checks, 4/5 completeness, 432/2,000 triple checks, exact **degree-7** linear ideal witnesses (`4` in $d=6$, `18` in $d=10$) & **degree-10** reduced SOS norm | ~18s |
| `04` | `scripts/04_clifford_6_8_12_mub_QED.mjs` | **Tier 2 (Clifford $\mathrm{Cl}_{q+1} \otimes \mathrm{MUB}(\mathbb{F}_q)$):** Breaks quaternionic $m \le 5$ ceiling — exact **6-MUB in $d=10$** ($M_{40}$, rank 4), **8-MUB in $d=14$** (rank 8), **12-MUB in $d=22$** (rank 32) | ~72s |
| `05` | `scripts/05_universal_jordan_wigner_N_Cstar_inf_QED.mjs` | **Tier 3 (Universal Jordan–Wigner Theorem, $N_{C^*}(d) = \infty$):** Exact finite-dimensional $m$-MUBs for every $d \ge 2$ and arbitrary $m > d+1$, verified over $\mathbb{Z}[\zeta_{2d}]$ for `(3,5)` ($9\times 9$), `(3,7)`, `(6,5)` ($36\times 36$), `(6,8)`, `(10,12)`, `(10,21)` | ~0.6s |
| `06` | `scripts/06_B_new_no_4th_vec_nullstellensatz_QED.mjs` | **Paper I:** $B_{\mathrm{new}} \in \frac{1}{6}\mathbb{Z}[\zeta_{20}][\sqrt{6}]$ (36 bases) Degree-6 Hilbert Nullstellensatz full rank `1256/1256` mod $p_1 = 1000081, p_2 = 1000921$ | ~6s |
| `07` | `scripts/07_B_new_fixed_Hconf_isolation_QED.mjs` | **Paper I:** Strict fixed-$\{I, H_{\mathrm{conf}}\}$ isolation of $B_{\mathrm{new}}$ (`90/90` full rank) & strict `4/4` transversality of $\mathcal{M}_4$ | ~0.1s |
| `08` | `scripts/08_B_iso_72orbits_171D_rigidity_QED.mjs` | **Paper I:** Chiral $B_{\mathrm{iso}} \in \frac{1}{9}\mathbb{Z}[i,\sqrt{5}]$ (`{0,0,16}`) 72-basis orbit & exact `171/171` Joint 3-MUB Rigidity | ~4s |

All file checksums are recorded in `PACKAGE_SHA256.json`.

---

## Citation

If you use or verify these results, please cite the archival release:

```bibtex
@misc{rui2026cstarmub,
  author       = {Rui, Yinghao},
  title        = {{Higher-Rank C*-Algebra MUBs ($N_{C^*}(d) = \infty$), Gribling--Polak Def. 3.1 Alignment, and $d=10$ Zero-Float Verification Suite}},
  month        = nov,
  year         = 2026,
  publisher    = {Zenodo},
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.24530887},
  url          = {https://doi.org/10.5281/zenodo.24530887}
}
```
