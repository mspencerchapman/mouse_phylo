# Figure files

Written when the Rmds are knitted with `resave_plots=TRUE` (the default; set it to FALSE in the setup chunk to knit
without overwriting). The two per-mouse Rmds take a `sample_ID` parameter, e.g.
`rmarkdown::render("Mouse_phylo_clustering.Rmd", params=list(sample_ID="MD7635"), output_file="Mouse_phylo_clustering_MD7635.html")`.
Panel numbers follow the manuscript draft of 2026-10-31. MD7634 = mouse 1, MD7635 = mouse 2.

| Panel | File(s) | Made by |
|---|---|---|
| Fig 1B-C, corrected ultrametric trees | `<mouse>_tree_corrected_ultrametric.pdf` | `Mouse_phylo_clustering.Rmd` |
| Fig 1D-E, chi-squared p-value by molecular-time cut-off | `<mouse>_chisq_by_moltime_plot.pdf` | `Mouse_phylo_clustering.Rmd` |
| Fig S1A-B, corrected trees (not ultrametric) | `<mouse>_tree_corrected.pdf` | `Mouse_phylo_clustering.Rmd` |
| Fig S1C-D, AMOVA on the corrected trees | `<mouse>_AMOVA_corrected_tree.pdf` | `Mouse_phylo_clustering.Rmd` / `Mouse_phylo_initial.Rmd` |
| Fig S1E-F, private-mutation VAFs and clonal / sub-clonal assignment | `<mouse>_private_VAF_by_depth.pdf`, `<mouse>_private_P_clonal.pdf` | `Mouse_phylo_initial.Rmd` |
| Fig S1E-F, signatures of shared, clonal and sub-clonal mutations | `<mouse>_signatures_by_clonality.pdf`; both mice pooled: `../burden_correction/results/signatures/signatures_by_clonality.pdf` | `Mouse_phylo_initial.Rmd`; `burden_correction/signatures_by_clonality.R` |
| Fig S[X], migration model | `migration_model/*.pdf` (below) | `migration_model/Migration_model.Rmd` |
| Not for the manuscript | `<mouse>_AMOVA_uncorrected_tree.pdf` (pipeline tree, for comparison); `amova/` (standalone four-panel AMOVA run) | |

Migration model (`migration_model/`):

- `trees_prior_range_T.pdf`: both trees with the 20-70 prior range of T shaded.
- `heatmap_<MD7634|MD7635|joint>.pdf` (panels A-B): posterior over T x migration rate (joint and given T), with the
  marginal posteriors of T (below) and of the rate (beside the joint posterior, with P(no migration)); colonisation
  windows marked.
- `ridges_<MD7634|MD7635|joint>.pdf`: migration rate given T, slices at T = 20, 25, ..., 70.
- `marginal_migration_rate.pdf` (panel E): posterior of the migration rate for mouse 1, mouse 2 and the joint fit,
  averaged over T, with the posterior probability of no migration as a bar at 0 and in the panel titles.
- `relative_likelihood_colonisation_window.pdf`: relative likelihood of the rate for T in the colonisation window (no
  longer in the figure; the text no longer conditions on the colonisation window).
- `posterior_predictive_chisq_curve.pdf` (panel F): observed chi-squared curve (all lineages, single colonies
  included) vs simulations from the posterior and with no clustering.
- The synthetic-likelihood (chi-squared only) fits were dropped on 2026-11-03; old PDFs in
  `../archive/synthetic_likelihood_plots/`.
- `posterior_predictive_same_bone_pairs.pdf`: observed vs simulated same-bone fraction of colony pairs.
