# Figure 1 — behavior & response sparsening

| Panel | Content | Script |
|-------|---------|--------|
| a, b | Setup, task and experimental-design schematics | — (Illustrator; `_schematic/`) |
| c | Training time course, example mouse (F07): Hit/FA/d′ | `behavior/code/plotBehaviorSummary_F07.m` |
| d | Trained-cohort performance in Hit-vs-FA (d′) space, D1–D6 | `behavior/code/plotBehaviorSummary_F07.m` |
| e | Example-session rasters (neurons × time), Naive vs Trained, Go/NoGo | `sparsening/code/sparsening_popmean_make.m` |
| f | Example population-mean timecourses (same sessions as e) | `sparsening/code/sparsening_popmean_make.m` |
| g | Grand-average population mean with significance bars | `sparsening/code/sparsening_popmean_make.m` |
| h | Fraction of responsive neurons per 1 s window, Naive vs Trained | `sparsening/code/fractionActive_make.m` |
| S14 | σ-floor diagnostic | `Supp_sigma_floor/code/make_sigma_floor_figure.m` |

Plot-data CSVs and stats: `behavior/code/behavior_make.m` plus the make scripts above.
Cache builders: `behavior/code/makeF07TrainingCSV.m` (from the public deposit),
`sparsening/code/sparsening_popmean_pull.m`, `sparsening/code/fractionActive_pull.m`.
