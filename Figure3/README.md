# Figure 3 — native-space geometry & dimensionality

| Panel | Content | Script |
|-------|---------|--------|
| a | Peak centroid distance by Go/NoGo angle (box+swarm) | `nativeSpaceGeometry/code/make_originalSpaceGeometryPlotting.m` |
| b–d | Centroid-distance timecourses: per angle, group difference, pooled + onset | `nativeSpaceGeometry/code/make_originalSpaceGeometryPlotting.m` |
| e–h | Same as a–d for the point-cloud radius | `nativeSpaceGeometry/code/make_originalSpaceGeometryPlotting.m` |
| dimensionality panels | Eigenspectrum with D_PR/k80/tail-slope overlays; D_PR timecourse (example angle + grand average with onsets); Trained−Naive difference for all three measures | `dimensionality/code/run_dim_production.m` (renders `make_plotDimSVMSummaries.m`) |
| S4 | Per-day native-geometry timecourse grids (all measures) | `nativeSpaceGeometry/code/make_perDay_combined_sheet.m` |
| S5 | Per-day dimensionality grids (D_PR, k80, tail slope) | `dimensionality/code/make_perDay_dimensionality_sheet.m` |

Stats and plot-data CSV export: `nativeSpaceGeometry/code/export_native_stats_csv.m`,
`dimensionality/code/make_dimensionality_windowed_stats.m` (includes the pooled-classes
"overall" D_PR reported in the text, from `dimensionality/code/pull_pooledSpectralPR.m`).
Cache builders: `nativeSpaceGeometry/code/pull_originalSpaceBasicGeometry.m`,
`dimensionality/code/pull_SVMPCATimeCourse_updated.m`, `dimensionality/code/pull_pooledSpectralPR.m`.
