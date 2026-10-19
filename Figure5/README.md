# Figure 5 — behavioral-state controls on capacity

| Panel | Content | Script |
|-------|---------|--------|
| a | Example sessions: capacity and licking by lick tertile (T1/T3) | `behavioralControls/code/make_fig5_capacityLicking.m` |
| b | Pooled tertile capacity + licking at 45°, with matched naive level | `behavioralControls/code/make_fig5_capacityLicking.m` |
| c | ΔCapacity vs ΔLicking (T3−T1, 1–2 s), one dot per session | `behavioralControls/code/make_fig5_capacityLicking.m` |
| d | T3−T1 % change per split (lick/loco/pupil) vs the Trained−Naive reference | `behavioralControls/code/make_fig5_pctBoxes_combined.m` |
| S7 | Super-mouse capacity + reconstructed lick (All/Lick/No-lick) | `behavioralControls/code/make_fig5_supp_supermouse.m` |
| S8–S9 | Pupil splits (pre-stimulus and evoked windows) | `behavioralControls/code/make_fig5_supp_pupil.m` |
| S10 | Locomotion splits (T1/T3/equalized) + all-trials capacity | `behavioralControls/code/make_fig5_supp_loco.m` |

Inputs are the caches under `behavioralControls/cache/` (lick traces, loco/pupil QC and
splits) and the capacity caches in `../Figure4/manifoldCapacity/cache/capacityTimecourse/`.
