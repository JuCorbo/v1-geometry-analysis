# Figure 6 — geometry ↔ behavior

| Panel | Content | Script |
|-------|---------|--------|
| a | Example-session licking decision proxies (P(lick) + raster) | `timingDecision/code/plotMainCapacity.m` |
| b | Capacity divergence onset (bootstrap CI) vs the decision proxies by day | `timingDecision/code/plotMainCapacity.m` |
| c–e | Capacity vs d′: correlation timecourse, 1–2 s scatter, difficulty regressed out | `geometryVsDprime/code/buildCapacityDprimeMain.m` |
| f | Example sessions: NoGo trials on the Go→NoGo centroid axis (PaCMAP) | `noGoAxisPosition/code/exampleAxisPosition.m` |
| g–i | NoGo P(lick) vs axis position: examples, aggregated ranks, slope distributions vs shuffled null | `noGoAxisPosition/code/plotAxisPositionVsLick.m` |
| S11 a–b | Detection vs decision lickers (got50 examples + distribution) | `timingDecision/code/plotLickerTypes.m` |
| S11 c–e | Per-angle capacity-difference dynamics, their correlation (pooling justification), onsets of the other geometric measures vs proxies | `timingDecision/code/plotSuppMeasuresPooling.m` |
| S12 | PaCMAP spread-normalized separation vs d′ (same as c–e) | `geometryVsDprime/code/buildPacmapSepVsDprime.m` |
| S13 a, e | Raw axis positions per angle (PaCMAP and native) | `noGoAxisPosition/code/plotAxisPositionVsLick.m`, `nativeAxisDistance/code/plotNativeAxisVsLick.m` |
| S13 b–d | Native-space axis-position versions of Fig 6g–i | `nativeAxisDistance/code/plotNativeAxisVsLick.m` (example: `exampleNativeAxis.m`) |

Onset CIs (panel b, S11e): `timingDecision/code/dimensionalityOnset.m` →
`bootstrapOnsets.m` → `capacityBootstrap.m` → `data/onsetCI.mat`; reported numbers in
`stats/timing_stats.txt` via `writeTimingStats.m`.
Cache builders: `noGoAxisPosition/code/computeAxisPositionVsLick.m`,
`nativeAxisDistance/code/computeNativeAxisVsLick.m`,
`timingDecision/code/recomputeDivergence_fromCache.m` (lick-divergence upper bound from
the cached lick traces).
