---
editor_options: 
  markdown: 
    wrap: 72
---

# What is still missing, updated 2026-07-28

Audited the three current files against every promise made in the response
to reviewers. The prose is in good shape. **The gap is almost entirely
figures and tables: we promise four new or revised exhibits to the
reviewers and none of them are in the documents yet.**

Confirmed by file inspection: the five images embedded in the manuscript
are byte-identical to the 2026-07-21 original, so no figure has been
replaced or added.

**Every exhibit named below already exists in this repo and is current.**
The repo path is given for each one so nothing has to be hunted for or
regenerated. If you want to rebuild any of them, the generating script is
named too; all of them read saved `.rds` files and run in seconds without
touching `precompute.R`.

------------------------------------------------------------------------

## A. Promised to the reviewers, not in the documents

| \# | Promise in the response | Insert this file | Rebuild with |
|----|----|----|----|
| A1 | "Sensitivity figure. Added as Figure [N]" | `figures/tdm_weight_sensitivity.png` | `analysis/tdm_weight_sensitivity.R` |
| A2 | "Side-by-side model table. Added as Table [N]" | Table text is inside the response letter; the underlying survival numbers are `output/figure3_etf_survival_by_temp.csv` | `analysis/figure3_tdm_curves.R` |
| A3 | Figure 5 caption promises numeric labels above each bar | `figures/mcda_composite_scores.png` — **already has the labels, and already reflects the current run.** Swapping it in fixes the stale bar values at the same time | `analysis/mcda.R` |
| A4 | "we now reproduce the individual allocations and the written justifications in full" | `output/table_S2-7_panelist_weights.csv` (Table S2-7) and `output/table_S2-8_justifications.md` (Table S2-8) | `analysis/elicitation_tables.R` |

Supporting numbers, if a caption needs them: `output/tdm_weight_sensitivity.csv`
and `output/tdm_weight_sensitivity_martin_sweep.csv` (A1),
`output/mcda_composite_scores_v2.csv` (A3).

## B. Internal inconsistencies that a reviewer would catch

| \# | Issue | Where the fix lives |
|----|----|----|
| B1 | **Figure 3 was not replotted.** The caption says 12.14 °C but the image is still the original single-panel daily-hazard plot — no 10–18 °C range, no cumulative-survival panel, no shaded operational band. Panel (b) is what actually rebuts "TDM.1 has essentially no temperature response". | `figures/figure3_tdm_curves.png`, from `analysis/figure3_tdm_curves.R` |
| B2 | **SI table numbering collides.** The parameter table was renumbered to S2-9, but Section S2.8 still says "presented in Table S2-7", and S2-7 is now the panelist table. Two `[EDITORIAL]` notes flag this and need removing once fixed. | Edit in the SI directly |
| B3 | **Anderson et al. (2022) and Bloomer et al. (2023) are cited but missing from the reference list.** Martin et al. (2020) is also needed once the verbatim justifications go in. | Formatted entries in `output/new_references.md` |
| B4 | Response says the panel weighting moved "within Section 2.5.1". It is actually a new **Section 2.5.2**. | Edit in the response |
| B5 | Response claims components not bearing on temperature were "condensed". Only one duplicated paragraph in 3.1.2 was struck. Either condense or soften the claim. | Edit in the response |
| B6 | **SI S2.8 claims the model "captured the magnitude and general temporal dynamics… including major peaks in 2015, 2018, and 2023–2024".** It does not — see D1. 2023 and 2024 are the two worst-fit years in the record, under by 20,152 and 28,599 spawners. | `output/calibration_predictions.csv` |

## C. Placeholders still in the text

| Where | What |
|----|----|
| Response | `[N]` ×3 and `[S_]` ×2 — resolve once the figures and tables are in |
| Response | `[Consequence Models]` — remove the brackets, the rename is done |
| Response | All `[INTERNAL: ...]` notes — delete before submission |
| SI | `[citation to be added]` for the Hazel Avenue statistics → the WTMP deliverable |
| SI | `[EDITORIAL: ...]` note on Table S2-9 numbering — delete once B2 is resolved |
| SI | `[EDITORIAL — DECISION NEEDED]` note on the fit statistics — replace per D1 |
| Manuscript | `[*add repository URL upon acceptance*]` — normal, leave until acceptance |

## D. Decisions that need a human

### D1. Calibration fit statistics — **RESOLVED, and the answer is bad news**

*Script:* `analysis/calibration_fit_statistics.R`
*Outputs:* `output/calibration_fit_statistics.csv`, `output/calibration_predictions.csv`,
`figures/calibration_observed_vs_predicted.png`, and `app_data/calib_pred_by_variant.rds`
(no longer 0 bytes)

The calibration-prediction step is restored and the statistics are now
reproducible from saved artifacts, without re-running `precompute.R`.

**The reproduction is exact.** Re-running the optimiser on the reconstructed
objective returns the saved parameters to machine precision — SAR 0.0026852 and
`rear_surv` 0.5427946, difference 0.00e+00 on both. So the numbers below are
what this model actually does.

| Statistic, TDM-weighted, 2014–2024 | SI currently says | Regenerated |
|---|---|---|
| R² (squared correlation) | 0.72 | **0.13** |
| R² (Nash–Sutcliffe) | — | **−0.59** |
| RMSE | 8,240 | **13,345** |
| MAPE | 24% | **51%** |

**The published values cannot be recovered under any variant tested.** Not with
all 14 years instead of 11 (R² = 0.31, RMSE = 11,829). Not per-variant. And not
by reverting the carrying capacity to the pre-rerun 12,493, which gives R² = 0.064,
RMSE = 11,812, MAPE = 53%. They belong to a structurally different earlier version
of the pipeline, not to a different parameterisation of this one.

A negative Nash–Sutcliffe means the model predicts observed escapement **worse
than the mean of the observations would**. The year-by-year residuals show why:

| Year | Observed | Predicted | Residual |
|---|---|---|---|
| 2016 | 14,473 | 34,148 | +19,675 |
| 2023 | 37,321 | 17,169 | −20,152 |
| 2024 | 45,541 | 16,942 | **−28,599** |

The model produces a nearly flat 17,000–25,000 trajectory and misses the
2023–2024 surge entirely.

**What to do — pick one:**

1. **Replace the sentence with the regenerated statistics and reframe the claim.**
   Honest and defensible if paired with the framing we already adopted for the
   projections: this is a *relative* comparison across alternatives under a fixed
   climatology, not a predictive hindcast. Two parameters are fitted to reproduce
   general magnitude; interannual variation is driven by ocean conditions the
   model does not represent. **This is my recommendation.**
2. **Delete the fit statistics and the "captured the temporal dynamics" claim
   (B6) entirely**, and describe the calibration as a magnitude-matching exercise.
   Safest, loses the least defensible content.
3. Do not quote R² = 0.72. It is not reproducible and the code is public.

### D2–D7, unchanged

| \# | Item |
|----|----|
| D2 | **Panelist 5's weights**: prose says 45/20/35, recorded row says 40/30/30. Visible the moment Table S2-8 reproduces the justification verbatim. Source: `output/table_S2-8_justifications.md` |
| D3 | **The 0.05%–0.68% SAR range** was replaced with verified brood-year values. Not reproducible from `app_data/SAR LAR Releases.xlsx` under any of twelve aggregations. If it came from another dataset, cite that. |
| D4 | **Acknowledgments** omit one panel member and spell another differently from the scoresheet. |
| D5 | **Anonymise or name the panelists** in Tables S2-7/S2-8. If named, we need consent; the justifications are candid. |
| D6 | **Three highlight colours** now layer in the manuscript (yellow, cyan, green) with no key. Green marks the 2026-07-28 automated edits. Collapse or add a legend before this goes to Brian or the editor. |
| D7 | **Title scoping.** Reviewer 1 objects that "Balancing" implies absolute performance. Conclusions were scoped; the title was left alone deliberately. |

------------------------------------------------------------------------

## Where everything lives

| | Path |
|----|----|
| Figures for the manuscript and SI | `figures/` |
| Tables and supporting numbers | `output/` |
| Scripts that regenerate them | `analysis/` |
| Full findings, organised as edits | `output/MANUSCRIPT_REVISION_HANDOFF.md` |
| Full findings, as an analysis record | `output/REVISION_FINDINGS.md` |
| Options memo for the app data question | `output/APP_DATA_UPDATE_OPTIONS.md` |
| Pre-edit copies of the three documents | `manuscript_backup/` |

**Note:** the live `.docx` files are no longer in the repo root — only the
pre-edit backups are. `.gitignore` excludes `*.docx`, so the working copies live
alongside the repo by design.

## Suggested order

1.  **A4 and A2 first.** Both are tables and both are copy-paste from files that
    already exist. A4 is the direct answer to the elicitation-rigour complaint.
2.  **B1, A1, A3.** All three images are rendered and current in `figures/`;
    they need inserting and captioning. A3 also fixes the stale bar values.
3.  **D1.** Decide before anything else goes to the editor — this is now a
    known-unreproducible number in a manuscript whose code is public.
4.  **B6, B3, B2, B4, C.** Cleanup, with B6 falling out of whatever D1 decides.

Once A1 through A4 are in, every `[N]` and `[S_]` in the response can be filled
with a real number, and the response stops making claims the manuscript does not
support.
