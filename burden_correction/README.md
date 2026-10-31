# Branch-length correction for sensitivity and sub-clonal mutations

Replaces the clonal/sub-clonal correction in `Mouse_phylo_clustering.Rmd` (old chunks `branch_level_clonality_analysis`
and `private_branch_correction`; now chunk `burden_correction`, since 2026-10-29). Developed 2026-10-27. The previous
Rmd, its chi-squared plots and figures are in `../archive/*_before_burden_correction*`.
`Mouse_phylo_initial.Rmd` was updated the same way on 2026-10-30. Its old binomial-mixture chunks (`binom_mix_model`,
`branch_level_clonality_analysis`, `private_branch_correction`) are replaced by `burden_correction`, `signatures_by_clonality`
and `private_branch_correction`. The previous version is `../archive/Mouse_phylo_initial_before_burden_correction.Rmd` / `.html`.

## Why the Rmd correction was replaced

- **Detection depends on depth.** Every private mutation on the trees (`*_postMS_reduced_a_j_vaf_post_mix_post_dup`) has
  NV >= 3 and VAF >= 0.2 (the tree-building `vaf` filter; 3 exceptions in ~26,000). Low-VAF (sub-clonal) mutations are
  removed far more often at low depth, so the sub-clonal fraction among called mutations falls with depth
  (clonal fraction 1.0 at <= 10x, ~0.47 at 18-25x). The Rmd took its prior from NR > 20 mutations (a minority; median
  sample depth is ~17-19x) and shrank every sample towards it.
- **The "0.33" sub-clonal peak** is a ~0.25-0.28 component truncated at VAF 0.2.
- **Male X/Y mutations** are hemizygous (VAF ~1, half depth): ~1/3 of all private mutations at <= 10x. They were treated
  as diploid.
- **Calls with exactly 3 variant reads are only partly detected** (~0.7 relative to NV >= 4; caller/LCM filters), which
  also makes low-depth mutations look clonal.
- The Rmd chunks also depended on objects only created in `Mouse_phylo_initial.Rmd`.

## Method (`R/burden_correction.R`)

SNVs only (the input trees were built with SNVs only, `-j`). For each individual:

1. Private SNVs of each sample: `get_private()`; X/Y outside the PAR get copy factor 2 (`copy_factor_mouse()`).
2. `fit_mixture()`: two-component **beta-binomial** mixture **conditioned on detection**. For a mutation with depth NR,

   P(NV | detected) = sum_k w_k BB(NV | NR, cf p_k, rho_k) d(NV) / sum_k w_k P(pass | NR, cf p_k, rho_k)

   with clonal VAF `pc`, sub-clonal VAF `ps` and overdispersion `rho_c`, `rho_s` shared across samples, the true-space
   clonal weight `w` per sample (profiled out), and detection weights `d(NV)` = `delta` for NV = 3, 4 (1 above).
   The clonal VAF is shared by all colonies (single-cell sorted; impure colonies should be removed beforehand - the
   pipeline already excluded colonies with peak VAF < 0.4).
3. Caller sensitivity: the pipeline's germline-SNP pick-up rate (CaVEMan + LCM filters;
   `/nfs/casm/team273jn/ms56/my_programs/Sensitivity_analysis_from_SNPs.R`). Its output file is gone, but the input trees
   were built with it, so it is recovered exactly as tip SNV count / tip edge length (`caller_sensitivity_from_tree()`).
4. `clonal_sensitivity()`: caller sensitivity x P(clonal SNV passes the tree-building rule), the latter averaged over the
   sample's depth at the mutation sites of other branches (not selected on this sample); X/Y weighted 0.033. `delta` is
   set to 1 here (low-NV caller losses are already in the germline-SNP rate). With `caller_sensitivity = NULL` the
   sensitivity comes from the depth model alone, including `delta` (used in the simulation).
5. Corrected private branch = sum of per-mutation clonal probabilities / clonal sensitivity (`correct_private()`).
6. Shared branches rebuilt from raw SNV counts and divided by 1 - prod(1 - sensitivity) over their descendant samples
   (`corrected_tree()`).

`correct_branch_lengths(tree, COMB_mats, sex)` runs everything and returns the per-sample table, the fit and the
corrected tree.

## Scripts

- `simulate_correction.R [n_reps]`: simulation with known truth at each sample's real depths, four scenarios
  (matched model; spread of sub-clonal VAFs 0.125/0.25/0.375; beta-binomial overdispersion rho = 0.02; 10 low-VAF
  artefacts per colony), comparing the new method with the Rmd hard calls, sensitivity-only correction
  and raw counts. Output: `results/simulation/`.
- `run_correction.R`: applies the correction to MD7634 and MD7635. Output in `results/`: `corrected_burdens.csv`,
  `fitted_parameters.csv`, `clonal_weight_by_depth.csv`, `tree_<mouse>_corrected.tree`,
  `burden_correction_<mouse>.Rds`, `burden_correction_plots.pdf`; logs in `results/logs/`.

## Results (2026-10-27)

### Simulation (2 replicates x 4 scenarios x both mice' real depths; true clonal burden 35-80 per colony)

Bias of the mean corrected private burden / per-sample RMSE (% of truth):

| Scenario | New, beta-binomial | New, binomial | Rmd hard calls | Observed |
|---|---|---|---|---|
| Matched | -1.1 / 9.0 | -0.9 / 9.0 | -20.4 / 22.4 | +64 / 72 |
| Sub-clonal VAF spread | -1.3 / 9.6 | +5.9 / 10.7 | -8.0 / 12.1 | +50 / 57 |
| Overdispersed (rho 0.02) | +0.1 / 11.0 | -12.7 / 16.1 | -26.8 / 29.0 | +60 / 67 |
| Low-VAF artefacts | -2.1 / 8.7 | 0.0 / 8.1 | -21.0 / 23.1 | +67 / 74 |

The beta-binomial version is within +-2% in every scenario and its error does not depend on depth (Spearman -0.01 to
0.09; Rmd 0.15-0.27). Fitted parameters recover the truth (pc 0.50, ps ~0.27, delta ~0.72). Per-sample RMSE ~9-11% is
mostly sampling noise. The simulation uses the depth-model sensitivity (caller losses = delta only); on the real
data the caller sensitivity is taken from the pipeline instead, so that part is not tested here. An earlier test with colony-to-colony clonal VAF variation (impure colonies) broke the method;
that scenario was dropped because the colonies are single-cell derived and impure ones are filtered out.

### Real data

| | MD7634 | MD7635 |
|---|---|---|
| pc / ps | 0.497 / 0.265 | 0.505 / 0.259 |
| rho clonal / sub-clonal | 0.016 / 0.006 | 0.007 / 0.013 |
| delta NV=3 / NV=4 | 0.74 / 0.99 | 0.68 / 1.00 |
| Caller sensitivity (median, range) | 0.844 (0.75-0.87) | 0.921 (0.85-0.94) |
| Clonal sensitivity (median) | 0.82 | 0.91 |
| Clonal weight per colony (median, range) | 0.60 (0.09-1.0) | 0.48 (0.09-0.93) |
| Private SNVs: raw / pipeline tree / corrected | 84.3 / 100.3 / 61.9 | 90.0 / 98.2 / 53.3 |
| Total burden: pipeline tree / corrected (SEM) | 130.3 / 91.9 (1.5) | 124.7 / 79.8 (2.0) |
| Spearman corrected private vs depth | 0.09 | -0.05 |

By tissue (corrected total): MD7634 femur 95.7, iliac 85.8, spine 91.9, tibia 94.8 (private burden differs by tissue,
Kruskal-Wallis p = 0.011, driven by iliac); MD7635 femur 81.9, iliac 80.0, spine 83.1 (n = 2), tibia 77.8 (p = 0.2).

### Open issues

- The clonal weight refitted within depth bins still drifts (MD7634 0.60 at <= 10x to 0.47 at 17-20x; MD7635 0.53 to
  0.40), which the simulation does not show, so something depth-related is not modelled (e.g. caller detection that
  depends on depth beyond the NV = 3 effect). The effect on corrected burdens is limited (no correlation with depth).
- The two mice differ by ~12 mutations in corrected total burden, mostly from caller sensitivity (0.84 vs 0.92) and the
  sub-clonal fraction; the earlier Rmd numbers (uniform 69.9 / 66.6; sample-specific 76.5 / 70.4) are not comparable.

## Signatures of clonal vs sub-clonal private mutations (2026-10-30)

`signatures_by_clonality.R` assigns each private SNV using its posterior P(clonal) from the mixture: weighted, most
likely assignment (MAP), and high-confidence (P > 0.9 / P < 0.1). It then builds 96-channel profiles and compares
them with shared (in vivo) mutations. Contexts were fetched read-only from the farm reference
(`/nfs/cancer_ref02/Mus_musculus/GRCm38/genome.fa`) and cached in `results/signatures/contexts.Rds`. Output is in
`results/signatures/`, including per-mutation assignments in `private_mutation_assignments.csv`.

Both mice combined:

| Set | n | C>T at CpG | SBS1 / SBS5 / SBS18 (restricted fit) | Cosine to shared |
|---|---|---|---|---|
| Shared | 1,822 | 0.21 | 0.23 / 0.62 / 0.16 | 1 |
| Clonal, P > 0.9 | 7,714 | 0.24 | 0.26 / 0.63 / 0.12 | 0.97 |
| Clonal, weighted | 14,110 | 0.20 | 0.21 / 0.67 / 0.13 | 0.96 |
| Sub-clonal, weighted | 10,847 | 0.10 | 0.07 / 0.79 / 0.14 | 0.79 |
| Sub-clonal, P < 0.1 | 4,337 | 0.05 | 0.02 / 0.84 / 0.15 | 0.66 |

Clonal private mutations match shared mutations (SBS1-rich, in vivo). Sub-clonal ones are SBS1-depleted, with flat
non-CpG C>T (a strict COSMIC fit picks SBS30-like). This fits in-culture mutations or low-VAF artefacts. Separation
grows with assignment confidence, which supports the mixture classification.
