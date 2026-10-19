# Figure 2 — population geometry (PaCMAP) & cross-method validation

| Panel | Content | Script |
|-------|---------|--------|
| a | Example sessions (Trained/Naive, 30°): PaCMAP cloud snapshots at 4 times | `pacmap/code/s03_render_cohort_and_theta_D135.m` |
| b | Example centroid trajectories (same sessions), θ indicated | `pacmap/code/s03_render_cohort_and_theta_D135.m` |
| c | Cohort PaCMAP, odd days (90/30/20°): aligned trajectories + clouds | `pacmap/code/s03_render_cohort_and_theta_D135.m` |
| d | θ trajectory-angle summary by day, Naive vs Trained | `pacmap/code/s03_render_cohort_and_theta_D135.m` |
| e | Overlaid Go/NoGo 1-SD ellipses, D1→D6, Naive vs Trained | `pacmap/code/s04_render_overlay_ellipses.m` |
| S1 | Cohort PaCMAP, even days (45/25/15°) | `pacmap/code/s03b_render_cohort_D246.m` |
| S2 | Same panels under baseline z-score | `s03`/`s03b`/`s04` with `PAC_BASEROOT` = `pacmap/cache/baselineZ` |
| S3 | Clouds + ellipses across 5 DR methods (PaCMAP/UMATO/cMDS/Isomap/nMDS) | `multiDR/code/s01_multiDR_embed_baselinez.m` → `s02_multiDR_figure_baselinez.m` |

θ stats and CSV export: `pacmap/code/s05_theta_stats_export.m` (reads
`pacmap/data/theta_bwgz_s0.csv`, which ships with the package and is rewritten
by `s03`).

Normalization per panel set: main + S1 = global per-neuron z-score with per-trial baseline
subtraction; S2 and S3 = baseline z-score (`multiDR/code/s03_knn_shortcircuit_analysis.m`
quantifies the neighbour-graph behaviour under each; see Methods, Cross-method validation).

Embedding pipeline (cached; rerun needs Python `pacmap==0.8.0`, see `pacmap/code/ENVIRONMENT.txt`):
`pacmap/code/s01_pull_data4timeDR_globalz_sub.m` / `s01b_pull_data4timeDR_baselinez.m` →
`s02_make_pacmap_umato.py` → per-session embedding caches.
