# Analysis code & data — "Geometric reshaping of task-relevant representations in the primary visual cortex supports perceptual decisions"

Self-contained reproduction package: every figure panel, plot-data CSV, and statistic
in the paper can be regenerated from the code and cached intermediates in this folder,
with no database or network access. Each figure folder has its own README describing
its modules; the manuscript's Methods section describes the analyses.
First versions of every script were initially written by Leyla Caglar and Julien Corbo,
then refined and checked with the help of LLMs.

## System requirements

- **MATLAB** R2023a or later with the Statistics and Machine Learning Toolbox.
  Windows, macOS, or Linux; no non-standard hardware. Tested on Windows 11 / R2023a.
- Nothing else to **reproduce the figures/stats** — all intermediates are included.
- Only to **regenerate intermediates from raw data** (optional):
  - the public dataset deposit (Zenodo, https://zenodo.org/records/8109858 — this study
    uses its Trained and Naive datasets; the Audiovisual dataset is not used);
  - **IBM CPLEX** + the `dnn-object-manifolds` toolbox (Chung, Lee & Sompolinsky 2018)
    for manifold capacity (`_capacity_source/`);
  - **Python** with `pacmap==0.8.0` (pinned) for the PaCMAP embeddings (Figure 2).

## Layout

```
Figure<N>/
  <analysis>/            one folder per distinct analysis feeding the figure
    code/                make_* scripts: load cache/ -> write panels/, data/, stats/
    cache/               cached intermediates the make scripts run from (included)
    panels/              rendered figure parts (PDF/PNG), main/ and supp/
    data/                tidy CSVs of the values actually plotted
    stats/               tests, n, statistics, exact p, effect sizes (txt)
  _schematic/            hand-drawn panels (setup/design) with no code or data
_sharedWindowedData/     the ONE shared windowed intermediate (per-session *_win.mat)
                         + accessor getWindowed.m + groupSeparationTest.m (cluster-based
                         permutation test used for all time-resolved group comparisons)
_sharedBehavior/         slim trial tables extracted from the public deposit
_capacity_source/        manifold-capacity generators (need CPLEX) + the bundled toolbox
```

## Installation guide

Clone the repository or download it as a zip and extract (~1 GB; the download is
the only time-consuming step — there is nothing to install or compile beyond a
standard MATLAB installation, and no paths to configure). Typical time on a normal
desktop: a few minutes, dominated by the download.

Windows note: some cache paths are long. If `git clone` reports "Filename too
long", enable long paths (`git config --global core.longpaths true`) or clone
into a short directory such as `C:\repro\`.

## Quick start

This package is the analysis code for the whole paper, not a standalone tool, so the
natural "demo" is reproducing one figure. In MATLAB, run
`Figure4/manifoldCapacity/code/make_capacityFullTimecourse.m`: it loads the included
capacity caches, runs the cluster-based permutation tests, and regenerates Figure 4's
manifold-capacity panels (`panels/main/main_capacity_fullData_timecourse.pdf` and the
per-day box-swarm), the plot-data CSV, and the stats file. Expected console output
ends with `Pooled: onset = 0.696 s` (cluster P = 1e-4) and the list of saved panels.
Run time: well under a minute on a normal desktop computer. Every other figure
regenerates the same way from its own `make_*` script.

## Instructions for use

Run any `make_*` (or numbered `sNN_*`) script in a `code/` folder from MATLAB; paths are
resolved relative to the script, so the package can live anywhere. Outputs overwrite
`panels/`, `data/`, `stats/` in place. Scripts print which cache file they load.
Every quantitative result in the manuscript reproduces this way; each figure folder's
README maps its panels to the scripts that produce them. Running every
figure-producing script in the package takes about 20 minutes in total on a
normal desktop.

## What is included vs. rebuilt

- All `cache/` intermediates are **included**, so every figure is reproducible offline.
- A few **builder** scripts document how an included intermediate was created from
  sources that are not distributed here (they are not needed for reproduction):
  - `_sharedWindowedData/code/make_windowed_cache.m` — builds the windowed cache from a
    staged CSV extraction of the public deposit (set `SCR` to your extraction folder);
    `validate_windowed_cache.m` proves it reproduces the database ground truth exactly.
  - `Figure1/behavior/code/makeF07TrainingCSV.m` — stages training behavior from the
    public deposit's DataFiles (set `zdir` to your download).
  - `_capacity_source/` — regenerating capacity caches requires IBM CPLEX.
  - Figure 2 `s01_*`/`s02_*` embedding steps require Python `pacmap==0.8.0`; the computed
    embeddings are cached, so the render/stats steps (`s03*`, `s04*`, `s05*`) run as-is.

## Reading the outputs

Keys needed to read the file names, CSVs and plots (all analysis details are in
the manuscript's Methods):

- **Days / Go/NoGo angles:** Go = 45° (fixed); NoGo per day D1–D6 = 135/90/75/70/65/60°,
  i.e. angular separation 90/45/30/25/20/15°.
- **Groups / colors:** Naive = gold `[0.85 0.65 0.13]` (n = 42 sessions, `BatuhanADNaive`),
  Trained = blue `[0.20 0.40 0.80]` (n = 60 sessions, `BatuhanAngDist`).
- Each `data/*.csv` is self-describing (header row); each `stats/*.txt` states the
  test, n, statistic, exact p-value, and the effect summary.

## License

MIT (see `LICENSE`). The bundled `dnn-object-manifolds` toolbox
(`Figure4/manifoldCapacity/_capacity_source/`) is by Chung, Lee & Sompolinsky
and retains its own MIT license.

