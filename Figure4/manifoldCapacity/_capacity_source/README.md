# Capacity computation — source scripts + toolbox

The manifold-capacity **figure deliverables** (`Figure4/manifoldCapacity`, `Figure5`, `Figure6/geometryVsDprime`, timing) are SQL-free: they read cached `*_capacity_timecourse.mat` files. Those caches were produced by the scripts here, which are **NOT SQL-free** — they are kept for provenance/reproducibility of how capacity was computed.

## generator/
- `rollingCapacity.m` — main generator. Flow: pulls `TraceSpikes` from the DB → builds the balanced `tuning_function` array `[nCells x nTrials x 2]` (Go/NoGo) per time window → saves `*_tuning.mat` → calls `check_capacity` (toolbox) → reads `capacity_results` and stores **capacity = 2 / capacity_results** in `*_capacity_timecourse.mat` (fields `capacity_rolling` [backward window] and `capacity_expanding`). SQL-based.
- `rollingCapacityLickNoLick.m` — lick vs no-lick split (Fig S7).
- `rollingCapacityResampling.m` — trial-count-matched / super-mouse resampling (Fig 5 / S7 controls).
- `rollingCapacityLickQuantiles.m` — lick-tertile (T1/T3) capacity, metric `LickAUCz_1_2` (Fig 5a–d caches).
- `rollingCapacityRandomTertileSplits.m` — random-thirds control capacity (`*_randThirdAvg_*`; Fig 5b/d reference).
- `rollingCapacityLickNoLickResampled.m` — super-mouse lick/no-lick FIXEDTRIALS summaries (Fig S7 caches).
- `locoCapacityPrep.m` + `locoCapacityCompute.m` — locomotion T1/T3/EQ split capacity (Fig S10 caches).
- `pupilCapacitySplits.m` (+ `runPupilCapacityFull.m` driver) — pupil T1/T3/EQ split capacity, baseline and post windows (Figs S8–S9 caches).

## toolbox_dnn-object-manifolds/  (MIT, sompolinsky-lab 2019)
Public Cohen/Chung/Sompolinsky **dnn-object-manifolds** MATLAB toolbox = the replica-based mean-field-theoretic manifold-capacity code cited in Methods (Chung, Lee & Sompolinsky 2018). Entry point used by the generator: `library/check_capacity.m`; point-cloud analysis in `point_cloud_analysis/manifold_stable_analysis_corr.m`; `FOptM/` is the optimization dependency. LICENSE preserved.

## Two generators
- `generator/rollingCapacity.m` (here) — the **original, SQL-based** generator (pulls `TraceSpikes` from the DB). Kept for provenance.
- `../code/pull_rollingCapacity.m` — the **SQL-free** equivalent used for the deliverable: same math, but the input comes from `getWindowed('sub')` instead of the DB. It genpath-adds this toolbox automatically.

## The one real dependency: IBM CPLEX (a solver, not data)
`check_capacity` needs no external data — it only reads the `*_tuning.mat` the generator writes. But its linear-separability test (`library/check_linear_seperability_svm_cplexqp.m`) calls **IBM CPLEX** (`cplexqp`/`cplexoptimset`). With CPLEX installed and on the path, `pull_rollingCapacity.m` runs **anywhere, SQL-free, end-to-end** (verified locally up to the CPLEX call: it loads the session, builds the tuning function, and enters the replica-MFT loop before needing the solver). Without CPLEX it errors at the separability step. CPLEX is free for academics via the IBM Academic Initiative. The toolbox ships **only** CPLEX separability checkers (no `quadprog`/`fitcsvm` fallback), so a CPLEX-free port would require swapping that solver and re-validating the capacity values.

## Note on the retained caches
The deliverable **reuses** the 324 retained `*_capacity_timecourse.mat` caches (capacity is computed from raw window sums, so it is independent of the baseline-z σ floor). `../code/make_capacityFullTimecourse.m` reads those and needs neither the DB nor CPLEX.
