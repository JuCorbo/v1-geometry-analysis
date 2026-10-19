# Supplementary — baseline-z σ-floor diagnostic

Whole-database diagnostic for the baseline-z sigma floor (σ_cap = 0.12; 34,973 neurons,
102 sessions). `code/make_sigma_floor_figure.m` reads the shared cache via
`getWindowed('baseline_nofloor')` and writes `panels/sigma_floor_diagnostic.{png,pdf}`.

Panels: (a) pooled σ distribution per group; (b) σ vs single-trial z SD (follows 1/σ below
σ ≈ 0.1, flattens above ≈ 0.12); (c) σ vs mean |z|; (d) per-session spread of the low-σ fraction.
