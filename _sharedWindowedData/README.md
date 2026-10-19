# Shared windowed data intermediate

The single pre-windowed intermediate that all neural analyses (Figs 1–4, and Fig 6 via
their caches) read, so that the windowing, baseline subtraction and σ rule live in exactly
one place.

## Contents

```
data/<group>/<day>/<DB>_<session>_win.mat     # group = trained | naive
manifest_trained.csv, manifest_naive.csv      # session, day, nTrials, nCells, nGo, nNoGo, path
code/getWindowed.m                            # THE accessor — use this, don't re-implement
code/groupSeparationTest.m                    # cluster-based permutation test (all time-resolved group comparisons)
code/make_windowed_cache.m                    # regenerator (from the public deposit; see below)
code/validate_windowed_cache.m                # ground-truth check
```

## How to read it

```matlab
S = load('.../BatuhanAngDist_F07_42221__win.mat');
[X, info] = getWindowed(S, 'baseline');     % X = [nTrials x nCells x 26], double
```

`scheme` is one of:

| scheme | what you get | used by |
|---|---|---|
| `raw` | raw integer window sums, no normalisation | Fig 1 fraction-active; Fig 4 capacity |
| `sub` | per-trial baseline subtracted (`X − 3b`), no z-score | — |
| `baseline` | baseline z-score, σ floored at 0.12 | Fig 1 pop-mean; Fig 2 supp; Fig 3; Fig 4 SVM |
| `global` | baseline-subtracted then global z-score | Fig 2 main |
| `baseline_nofloor` | as `baseline` with unfloored σ | σ-floor sensitivity analysis only |

`info` carries `tvec`, `tvec_centre`, `bins`, `centers`, `baselineBins`, trial metadata
(`trialID`, `Visual_Stim`, `Response`, `Outcome`), cell metadata (`cellID`, `OSI`, `PrefOri`),
and `provenance`.

## Conventions

- **Windows:** backward, non-overlapping, 3 frames, labelled by last frame, spanning −1 to 4 s:
  `bins = (c-2):c` for `c = 9:3:84` → `{7,8,9} … {82,83,84}`, 26 bins. `dt = 0.0646`,
  `rasterLen = 140`, stimulus onset at frame 25. Every bin is a full 3 frames.
- **Baseline bins** are 1–6 (frames 7–24); post-stimulus bins start exactly at frame 25.
- **Two time vectors:** `tvec` = last-frame time; `tvec_centre` = true window centre.
  Onsets are reported at `tvec_centre`.
- **σ floor = 0.12** on baseline z-scoring (rationale in `getWindowed.m` and the
  `Figure1/Supp_sigma_floor/` diagnostic).
- **Baseline subtraction** is recoverable because all bins are equal width:
  `subtracted = X − 3b`, with `b` stored per (trial, cell) in the file.
- The 1 s windows used by fraction-active are exact unions of 5 consecutive bins
  (baseline = bins 2–6; 0–1 s = 7–11; 1–2 s = 12–16; 2–3 s = 17–21).

## Validation

Fraction-active recomputed from this cache matches the values pulled from the original
SQL database (`Figure1/sparsening/cache/fractionActive_pulled.mat`) **exactly**:
trained 60/60 sessions and naive 42/42 sessions, max diff 0.000e+00. That single metric
validates decode, frame alignment, QC inclusion, trial labelling and binning together.
QC counts also match the source CellTables exactly (e.g. `F07_42221_` → 261 cells,
`056_20623_0` → 427 cells). `validate_windowed_cache.m` reruns the check.

## Limits

- **No frame-level traces.** Anything needing a different windowing must be regenerated
  from the deposit CSVs, not from this cache.
- Inclusion is fixed at the standard QC clause and Visual-block-only trials; changing QC
  means regenerating.
- Trials/cells are ordered by ascending ID; `cellID` and `trialID` are stored, so join back
  to the deposit tables by ID, not by position.

## Regenerating (not needed for reproduction)

```matlab
GROUP='trained';  make_windowed_cache      % ~11 min (after the CSV staging pass)
GROUP='naive';    make_windowed_cache
GROUP='trained';  validate_windowed_cache
```

Source: the public Zenodo deposit archive (`DatasetsCorbo2024.rar`) — `Trained_*` and
`Naive_*` tables (`AudioVisual_*` is a different cohort, not used). Members are individually
seekable with bsdtar (`tar.exe -xOf <rar> <member>`). The heavy `*_CellTrialTable.csv`
members (171/163 GB) are first streamed through the awk QC filter `code/stage1gen.awk`
(`-v nskip=1` trained, `-v nskip=2` naive — the naive export has an extra `layer` column),
reducing them to ~14/17 GB; budget ~95 min per group for that pass. The two exports have
different column layouts; `make_windowed_cache` handles this per group.
