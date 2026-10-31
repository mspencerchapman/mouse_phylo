# Localisation + migration model for anatomical clustering of HSC colonies

Question: the trees cluster by bone, which the paper reads as HSCs staying "predominantly in the same place". How often
do HSCs move between bones after they first localise? Developed 2026-10-30.

## Model (`R/migration_model.R`)

Time = molecular time on the corrected ultrametric tree (as in `Mouse_phylo_clustering.Rmd`: `make.ultrametric.tree` on
`burden_correction/results/tree_<mouse>_corrected.tree`, rescaled to the mean corrected burden; root 0, tips D ~ 80-90).

1. **Localisation at time T.** Every lineage crossing T is placed in a bone k with probability pi_k (the observed
   fraction of colonies from that bone), independently of the other lineages.
2. **Migration after T.** Each lineage has migration events at a constant rate mu; the destination is drawn from pi
   (F81 model). Real moves (to a different bone) happen at rate mu (1 - sum pi^2) on average.
   Results are reported as **expected moves per lineage between localisation and sampling** = mu (1 - sum pi^2)(D - T).
   The probability that a sampled HSC never moved since localisation is exp(-moves).

pi is the stationary distribution, so the likelihood of the observed bones is exact: the product, over lineages crossing
T, of sum_k pi_k L_k (Felsenstein pruning below T; `loglik_grid()`). This is checked against brute-force enumeration,
normalisation and forward simulation on a toy tree, and no-clustering limits on a 300-tip tree (`tests/test_likelihood.R`). A grid over T (1 mutation steps) x mu (0 plus log-spaced 10^-3.5 - 10^0.5) gives
the posterior.

- Priors: mu has 50% mass on 0 ("no migration", model 1) and 50% log-uniform over mu > 0 (model 2). T is uniform over
  5-50 or 1-D molecular time, or over a window calibrated to calendar time (below).
- Ultrametric tree: iteratively reweighted means (`make.ultrametric.tree`), not tip extension. Private branch lengths
  have to reflect real time (decision 2026-10-30). All tips are then at the mean corrected burden = age 6 months.
- Tree uncertainty: the likelihood is averaged over 100 parametric bootstrap trees (Poisson-resampled corrected branch
  lengths, made ultrametric again).
- Joint fit: the two mice share T as a fraction of tree depth and share mu.

## Rmd

`Migration_model.Rmd` (knits in ~20 s) is the main write-up. It uses a broad prior on T, uniform over 20-70
mutations, because the embryonic mutation rate is not constant and the molecular time of marrow colonisation is
uncertain. It shows, for each mouse and jointly:
- heatmaps of the posterior over T and the migration rate, both the joint density and the rate given T;
- P(no migration | T);
- posterior predictive checks and validation.

It also acknowledges that localisation is protracted, so the inferred migration rate is an upper bound on true
bone-to-bone migration. The calibrated-window results below came from `run_migration.R` and are superseded by the broad
prior in the Rmd.

Simulation-based section (added 2026-10-30). The model is fitted to the clade x bone chi-squared curve (cut-offs
10-50) by **synthetic likelihood** on a grid (`R/abc.R`; reference tables cached in `results/abc_reference_<mouse>.Rds`,
1,000 simulations per cell). Kernel ABC (rejection-style) is also included but was nearly uninformative, because a
single distance scale does not suit cells with very different noise. A relative-likelihood plot gives prior-free
statements about high migration rates.

## Script

`run_migration.R [n_boot_trees=100] [n_sim_reps=40]` (~5 min on the laptop). Output in `results/`:

- `posterior_summary.csv`, `posterior_summary_joint.csv`, `moves_given_T.csv`;
- `migration_model_plots.pdf` (likelihood surface; moves given a fixed T);
- `validation_simulations.csv` / `validation_summary.csv`;
- `posterior_predictive.csv` / `.pdf`;
- `migration_model_fits.Rds`;
- `logs/`.

## Results (2026-10-30)

### Validation

Bones were simulated on the real MD7634 tree with known T and moves, 40 replicates each, then fitted with the 5-50
prior.

| True T | True moves/lineage | Median T | 95% cover T | Median moves | 95% cover moves | Mean P(no migration) |
|---|---|---|---|---|---|---|
| 10 | 0 / 0.1 / 0.3 / 1 | 10 / 10 / 10 / 26 | 1.00 / 0.98 / 0.98 / 0.93 | 0 / 0.11 / 0.34 / 0.98 | 1.00 / 0.98 / 0.93 / 0.93 | 0.99 / 0 / 0 / 0 |
| 20 | 0 / 0.1 / 0.3 / 1 | 20 / 20 / 20 / 29 | 1.00 / 1.00 / 0.98 / 0.98 | 0 / 0.09 / 0.31 / 1.44 | 1.00 / 0.93 / 0.93 / 0.88 | 0.98 / 0 / 0 / 0 |
| 35 | 0 / 0.1 / 0.3 / 1 | 35 / 35 / 35 / 30 | 0.95 / 0.98 / 0.98 / 1.00 | 0 / 0.08 / 0.28 / 2.7 | 1.00 / 0.88 / 0.95 / 0.88 | 0.94 / 0.08 / 0.02 / 0 |

With up to ~0.3 moves per lineage, T and the move rate are recovered well. At ~1 move per lineage the clustering is
mostly erased, so T becomes poorly identified and the upper bound on moves is weak. With no migration and early
localisation the tree would be almost perfectly sorted by bone (posterior predictive "no_migration_T5"), which is far from
the data.

### Real data

The bone signal in the trees is weak in likelihood terms:

| | MD7634 | MD7635 |
|---|---|---|
| logL, bone independent of the tree | -213.1 | -153.5 |
| logL, best no-migration model | -213.1 (T = 74-80, i.e. no localisation) | -147.2 (T = 41-43) |
| logL, best model with migration | -210.0 | -141.8 |
| Fraction of same-bone pairs (random expectation) | 0.25 | 0.33 |
| Same-bone fraction, pairs with MRCA >= 20 / >= 30 | 0.34 (171 pairs) / 0.39 (85) | 0.41 (114) / 0.80 (39) |

Posterior (bootstrap-averaged trees, T prior 5-50):

| | MD7634 | MD7635 | Joint (point tree) |
|---|---|---|---|
| T (median, 95% CrI) | 27 (7-49) | 34 (27-44) | 31 (22-35) |
| Moves per lineage after T | 1.4 (0.40 - >100) | 0.10 (0 - 0.42) | 0.51 (0.30-0.98) |
| P(no migration) | 0 | 0.24 | 0 |
| P(never moved since T) | 0.24 (0 - 0.67) | 0.91 (0.66-1) | 0.60 (0.37-0.74) |

T and the move rate trade off. Fixing T gives these moves per lineage after T (median, 95% CrI):

| T | 5 | 10 | 15 | 20 | 25 | 30 | 40 |
|---|---|---|---|---|---|---|---|
| MD7634 | 14 (1.6 - >100) | 11 (1.4 - >100) | 3.2 (0.8 - >100) | 1.3 (0.6 - >100) | 1.0 (0.5 - >100) | 0.9 (0.46 - >100) | 2.7 (0.38 - >100) |
| MD7635 | 1.6 (0.8 - 100) | 1.3 (0.7 - 83) | 1.0 (0.5 - 48) | 0.56 (0.25 - 2.8) | 0.32 (0.15-0.73) | 0.17 (0.05-0.47) | 0 (0 - 0.10) |

Interpretation, conditional on the model:

- If localisation is early (T <~ 15-20 molecular time), both mice need roughly one or more moves per lineage
  afterwards. Lower 95% bounds are ~0.5-1.6, and the upper bounds are not informative. Otherwise the trees would be
  much more sorted by bone.
- MD7635 is best explained by later localisation (T ~ 30-35) followed by few moves (~0.1-0.4 per lineage). MD7634 has too
  little signal to separate "late localisation" from "early localisation with frequent migration". It shows only a
  modest excess of same-bone pairs at MRCA height ~18-35 and nothing later.
- Converting T and the move rate to calendar time needs a molecular-time to age mapping, in particular the molecular
  time at birth / marrow colonisation. The early (embryonic) mutation rate is higher, so this is not linear in age.
  See the calibration section below.

### Calendar calibration (Kapadia et al. 2025, Nature 641:681, doi:10.1038/s41586-025-08625-8)

C57BL/6J female HSC/MPP colonies:
- **before birth:** ~48.2 SNVs at birth (mixed-model intercept); ~10 by gastrulation (~E6.5-7). That gives ~2.9 SNVs/day
  before birth, assumed linear (rough).
- **after birth:** 45.3 SNVs/year (59.5 at 3 months), so a 6-month HSC has ~70.9 SNVs.

`mol_to_age()` / `age_to_mol()` implement this. Our corrected depths are higher (MD7634 91.3, MD7635 79.2; male mice,
different pipeline), so two versions are given:
- **unscaled:** our molecular time taken at face value;
- **scaled:** our tree depth mapped to 70.9 (scale 1.29 / 1.12).

Localisation window = bone-marrow colonisation, ~E15.5 to ~P7:
- unscaled: 36.6-49.1;
- scaled: MD7634 47.2-63.2, MD7635 40.9-54.9.

Posterior with the calibrated window (bootstrap-averaged trees). Moves are counted from localisation to sampling, which
is ~0.5 years:

| | MD7634 unscaled | MD7634 scaled | MD7635 unscaled | MD7635 scaled | Joint, scaled |
|---|---|---|---|---|---|
| P(no migration) | 0 | 0 | 0.74 | 0.80 | 0 |
| Moves per lineage | 2.5 (0.23 - >100) | 1.5 (0.04 - 72) | 0 (0 - 0.17) | 0 (0 - 0.11) | 0.21 (0.03 - 5.7) |
| Moves per year | 4.9 (0.47 - >100) | 2.9 (0.07 - >100) | 0 (0 - 0.34) | 0 (0 - 0.22) | 0.41 (0.06 - 11) |
| P(never moved since T) | 0.09 (0 - 0.79) | 0.23 (0 - 0.97) | 1 (0.84-1) | 1 (0.89-1) | 0.81 (0 - 0.97) |

Calendar ages of fixed T (from `moves_given_T.csv`): molecular time 20 ~ E9-10, 30 ~ E12.5-13.5, 40 ~ E15.5-16,
45 ~ E16.5-18.5, 50 ~ P15 unscaled / E18.5 scaled.

- **MD7635.** It is consistent with no migration after marrow colonisation (upper 95% bound ~0.2-0.3 moves per year).
  However, the likelihood peaks at T ~ 30 (~E12.5-13.5), 2.3 log-units above T = 40. So related lineages share a bone
  from before marrow colonisation, i.e. at fetal-liver stage. One reading is correlated seeding (sister lineages from
  the fetal liver seeding the same bone); the other is that the clock calibration is off by a few days.
- **MD7634.** It has almost no bone signal (+3 log-units over "no clustering"). With localisation at colonisation it
  needs several moves per year, but the interval is very wide.
- **The two mice disagree**, so the joint estimate (~0.4 moves per year, 0.06-11) mainly reflects that; report them
  separately.
- **Little late-tree information.** In the calibrated window only a few clades coalesce after T (MD7635: 39 pairs with
  MRCA >= 30, 9 with MRCA >= 40). Estimates of post-colonisation migration therefore rest on a handful of clades.

### Caveats

- **Migration model.** Migration is to any bone in proportion to pi, at a constant rate. Moves between neighbouring
  bones (femur-tibia, iliac-spine) or a burst of migration around localisation would look different. Any of these
  could explain the imperfect fit for MD7634, where the posterior predictive slightly under-predicts the same-bone
  excess at h ~ 18-35.
- **Sampling.** pi is the sampled fraction per bone (set by design), not the true HSC distribution; this matters
  only through the F81 destination probabilities. MD7635 has only 2 spine colonies.
- **Tree structure.** Lineages crossing T are treated as independent. Early polytomies and the ultrametric heights of
  shallow nodes (few mutations) are the main uncertainty in T; bootstrap averaging covers branch-length noise, not
  topology errors.
- **Interpretation of T.** "Not localised before T" means that lineages split before T are placed independently.
