# Figure 4 — separability: manifold capacity & linear decoding

| Panel | Content | Script |
|-------|---------|--------|
| a | Capacity by Go/NoGo angle (box+swarm, 0–3 s mean) | `manifoldCapacity/code/make_capacityFullTimecourse.m` |
| b–f | Capacity timecourses: example angle, per-angle (both groups), pooled + onset, per-angle difference, 1–3 s slope distribution | `manifoldCapacity/code/make_capacityFullTimecourse.m` |
| g–l | Same as a–f for linear-SVM decoding accuracy | `linearSeparability/code/make_svmFullTimecourse.m` |
| S6 | Per-angle capacity and SVM timecourse grids | same two make scripts (supp outputs) |

Cache builders: `linearSeparability/code/pull_svmAccuracy_fromCache.m` (SVM accuracy from
the shared windowed cache), `manifoldCapacity/code/pull_rollingCapacity.m` (regenerates the
per-session capacity caches in `manifoldCapacity/cache/capacityTimecourse/`; requires IBM
CPLEX — see `manifoldCapacity/_capacity_source/`).
