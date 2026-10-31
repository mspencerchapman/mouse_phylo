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
- `heatmap_<MD7634|MD7635|joint>.pdf`: posterior over T x migration rate (joint and given T), with the marginal posterior
  of T and P(no migration | T); colonisation windows marked.
- `ridges_<MD7634|MD7635|joint>.pdf`: migration rate given T, slices at T = 20, 25, ..., 70.
- `relative_likelihood_colonisation_window.pdf`: prior-free relative likelihood of the rate, T in the colonisation
  window (supports "excluded rates above ~0.2 moves per year" for mouse 2).
- `heatmap_*_synthetic_likelihood.pdf`, `ridges_*_synthetic_likelihood.pdf`: the same fits using only the chi-squared
  statistics ("fitting the model only to the chi-squared statistics gave consistent results").
- `posterior_predictive_chisq_curve.pdf`: observed chi-squared curve vs simulations from each posterior.
- `posterior_predictive_same_bone_pairs.pdf`: observed vs simulated same-bone fraction of colony pairs.
