---
editor_options: 
  markdown: 
    wrap: 72
---

# What is still missing, updated 2026-07-28

Audited the three current files against every promise made in the
response to reviewers. The prose is in good shape. **The gap is almost
entirely figures and tables: we promise four new or revised exhibits to
the reviewers and none of them are in the documents yet.**

Confirmed by file inspection: the five images embedded in the manuscript
are byte-identical to the 2026-07-21 original, so no figure has been
replaced or added.

**Every exhibit named below already exists in this repo and is
current.** The repo path is given for each one so nothing has to be
hunted for or regenerated. If you want to rebuild any of them, the
generating script is named too; all of them read saved `.rds` files and
run in seconds without touching `precompute.R`.

------------------------------------------------------------------------

## A. Promised to the reviewers, not in the documents

| \# | Promise in the response | Insert this file | Rebuild with |
|----|----|----|----|
| A1 | "Sensitivity figure. Added as Figure [N]" | `figures/tdm_weight_sensitivity.png` | `analysis/tdm_weight_sensitivity.R` |
| A2 | "Side-by-side model table. Added as Table [N]" | Table text is inside the response letter; the underlying survival numbers are `output/figure3_etf_survival_by_temp.csv` | `analysis/figure3_tdm_curves.R` |
| A3 | Figure 5 caption promises numeric labels above each bar | `figures/mcda_composite_scores.png` — **already has the labels, and already reflects the current run.** Swapping it in fixes the stale bar values at the same time | `analysis/mcda.R` |
| A4 | "we now reproduce the individual allocations and the written justifications in full" | `output/table_S2-7_panelist_weights.csv` (Table S2-7) and `output/table_S2-8_justifications.md` (Table S2-8) | `analysis/elicitation_tables.R` |

Supporting numbers, if a caption needs them:
`output/tdm_weight_sensitivity.csv` and
`output/tdm_weight_sensitivity_martin_sweep.csv` (A1),
`output/mcda_composite_scores_v2.csv` (A3).

------------------------------------------------------------------------

## A-bis. Where each exhibit goes, what number it gets, and what has to change

Numbering follows order of first mention, which is what the journal
copyeditor checks. Current first-mention order is Table 1 (§2.3), Table
2 (§2.4), Table 3 (§3.1), Figure 1 (§2.1), Figure 2 (§2.5), Figure 3
(§2.5.3), Figure 4 (§3.1.2), Figure 5 (§3.2).

| Exhibit | Number | Effect on existing numbering |
|----|----|----|
| TDM model comparison table (A2) | **Table 3** (new) | Current Table 3 becomes **Table 4** |
| TDM weight sensitivity figure (A1) | **Figure 6** (new) | None |
| Replotted TDM curves (B1) | Figure 3, unchanged | Caption becomes two-panel |
| MCDA composite scores (A3) | Figure 5, unchanged | None |
| Panelist weights (A4) | Table S2-7, already numbered | None |
| Panelist justifications (A4) | Table S2-8, already numbered | None |
| Front-loading decomposition (optional) | **Figure S2-1** (new) | None |
| Calibration observed vs predicted (D1) | **Figure S2-2** (new) | None |

### Table 3 (new) — TDM model comparison

Tables are numbered by first mention, and the comparison table is first
mentioned in §2.5.2, which precedes §3.1. It therefore takes the number
3, and the consequence table becomes Table 4.

**Insert:** at the end of new Section 2.5.2 (Temperature-Dependent
Mortality Model Weighting), immediately after the paragraph describing
the Delphi allocation and the final weights, so the table sits directly
above the weights it explains. That placement is what Reviewer 2 asked
for.

**Add this callout** at the end of that paragraph:

> Table 3 summarizes the three formulations by data source, run type,
> calibration range, temperature sensitivity, assigned weight and
> principal limitation.

**Caption:**

> **Table 3.** Comparison of the three temperature-dependent mortality
> (TDM) formulations evaluated for the lower American River, showing
> data source, functional form, run type represented, development
> context, calibration temperature range, relative stage sensitivity,
> principal limitation, and the weight assigned by expert elicitation
> (Section 2.5.2).

**Second callout, in the Discussion:** append "(Table 3)" to the
sentence identifying the 12.14 °C threshold as a winter-run field
estimate.

**Renumbering checklist — current Table 3 becomes Table 4.** Four
callouts in the manuscript and two in the SI, plus the caption itself:

| Document | Location | Current text |
|----|----|----|
| MS | §3.1, first paragraph | "…biological benefits and economic costs (Table 3)." |
| MS | §3.1.3 Hydropower Costs | "…scaled directly with bypass volume (Table 2, Table 3)…" |
| MS | §3.1.2 | "…worst to best alternative for Chinook salmon (Table 3)." |
| MS | §3.2 | "…little contribution from fish benefits (Table 3, Figure 5)." |
| SI | S2 | "…four meteorological scenarios as seen in Table 3 in the main text." |
| SI | S2 | "…in contrast to 7,600-11,073 seen in Table 3)." |

If these are Word cross-reference fields they update on F9; if plain
text, change all six by hand.

### Figure 6 (new) — TDM weight sensitivity

First mention falls in §3.2, after Figure 5, so no renumbering.

**Insert:** Section 3.2 (Trade-off Evaluation), immediately after the
paragraph reporting the composite scores and before the paragraph on
Reclamation's selection of PB4.

**Add this paragraph** as the callout:

> Because the TDM model weights were elicited rather than estimated, we
> evaluated how sensitive the composite ranking is to them (Figure 6).
> PB1 remained the top-ranked alternative under the Bratovich-only,
> Bartholow-only, equal and elicited weightings. Holding the other two
> formulations at their elicited ratio and increasing the weight on
> Martin et al. (2017) from 0 to 1, PB1 retained the top rank up to a
> weight of 0.988. The no-bypass alternative took the top rank only at a
> weight of 1.000, where the projected adult population index spans just
> 3.6 to 84.8 spawners across all nine alternatives, so the Chinook
> objective is normalized across a range of a few dozen fish and the
> ranking is decided almost entirely by the hydropower objective.

**Caption:**

> **Figure 6.** Sensitivity of the multi-criteria decision analysis to
> the expert-elicited temperature-dependent mortality model weights. (a)
> Composite score for each of the nine management alternatives under
> five weighting schemes: each formulation weighted alone, equal
> weighting, and the elicited weighting (0.51 Bratovich et al. 2020,
> 0.24 Bartholow and Heasley 2006, 0.25 Martin et al. 2017); the
> top-ranked alternative in each scheme is highlighted. (b) Composite
> score as the weight on Martin et al. (2017) is increased from 0 to 1
> with the other two formulations held at their elicited ratio.
> Objective weights are held at 0.40 Chinook salmon, 0.50 hydropower and
> 0.10 steelhead throughout, and meteorological scenarios at equal
> probability.

**Third callout, in the Discussion:** add "(Figure 6)" alongside the
existing robustness sentence.

### Figure 3 — replace image, rewrite caption

No number change. Swap the image and replace the caption, which
currently describes a single panel.

**Caption:**

> **Figure 3.** Temperature-dependent mortality (TDM) formulations for
> Chinook salmon early life stages. (a) Daily instantaneous mortality
> rate. The Martin et al. (2017) linear-threshold model (yellow) begins
> accruing mortality at 12.14 °C, the lowest threshold of the three
> formulations; the Bartholow and Heasley (2006) exponential model
> (green) and the Bratovich et al. (2020) exponential model (blue) begin
> at higher temperatures. Dashed lines indicate alevin-stage parameters.
> (b) Cumulative egg-to-fry survival at constant temperature, using
> accumulated thermal unit stage durations. The shaded band is the 5th
> to 95th percentile of modelled October to November daily water
> temperature at Hazel Avenue across all alternatives and meteorological
> years (14.1 to 18.2 °C). Dotted lines mark where each exponential
> formulation's cumulative survival falls below that of Martin et al.
> (2017), at 16.59 °C and 17.25 °C.

The existing callout in §2.5.3 needs no change, but add a pointer to
panel (b) in the sentence describing the three formulations, since panel
(b) is what carries the argument.

### Figure 5 — replace image only

No number or caption change. The caption already promises the numeric
bar labels, and the new image already has them and reflects the current
run.

### SI Tables S2-7 and S2-8

Both are already referenced in the SI text; only the tables are missing.

**Insert both** in Section S2.5, immediately after Table S2-6 and after
the paragraph beginning "Individual weightings are reported in Table
S2-7."

**Captions:**

> *Table S2-7. Individual expert weightings of the three
> temperature-dependent mortality formulations. Weights were allocated
> independently in a single round, each summing to 1.0. The panel mean
> is the weighting used in the analysis.*

> *Table S2-8. Written justifications recorded by each panelist
> alongside their weight allocation, reproduced verbatim.*

Resolve **D2 and D5 before inserting S2-8.** Panelist 5's prose states
45/20/35 against a recorded row of 40/30/30, which becomes visible the
moment the justification appears beside the table.

### SI cross-reference fix (B2)

Section S2.8 reads "Final calibrated parameter values are presented in
Table S2-7." That number now belongs to the panelist weighting table.
Change to **Table S2-9**, then delete the `[EDITORIAL]` note on the S2-9
caption.

### Figure S2-1 (optional) — front-loading decomposition

Worth adding because the Discussion now asserts the front-loading
mechanism as a result rather than a conjecture, and nothing in the paper
shows it.

**Insert:** Section S2.4, or a short new subsection following the TDM
formulations.

**Caption:**

> *Figure S2-1. Decomposition of the front-loading effect. (a) Daily
> water temperature at Hazel Avenue relative to the no-bypass
> alternative. (b) Egg-to-fry survival by spawn date under the Martin et
> al. (2017) formulation for the no-bypass alternative. (c) The same
> curves expressed as a difference from no bypass, showing November
> cohorts losing and October cohorts gaining. (d) Change in incubation
> mortality relative to no bypass, split at each alternative's crossover
> date into the cooling benefit earned before and the warming penalty
> paid after. (e) Cohort decomposition of the difference between PB6 and
> PB4.*

Panel (c) carries the argument: it shows the sign reversal by spawn date
and that December onward is exactly zero.

### Figure S2-2 — calibration observed versus predicted

Only under **D1 option 1**, that is, reporting the regenerated
statistics rather than deleting them.

**Insert:** Section S2.8, with the revised fit statistics.

**Caption:**

> *Figure S2-2. Observed versus predicted in-river escapement over the
> 2014–2024 calibration period, TDM-weighted. The model reproduces the
> general magnitude of escapement but not its interannual variation;
> residuals are largest in 2016, 2023 and 2024.*

Under D1 option 2, skip this figure. A plot showing a flat trajectory
against the 2023–2024 surge invites exactly the question that deleting
the statistics is meant to avoid.

### Filling the response placeholders

Once the above are in, replace throughout the response letter:

| Placeholder | Becomes |
|----|----|
| "Sensitivity figure. Added as Figure [N]" | Figure 6 |
| "Side-by-side model table. Added as Table [N]" | Table 3 |
| "new Table [N]" in the TDM characterization paragraph | Table 3 |
| "new Figure [N]" in the RSU response | Figure 6 |
| "Appendix [S\_]" (×2) | Section S2.5, Tables S2-7 and S2-8 |
| "[Consequence Models]" | Consequence Models, brackets removed |
| "within Section 2.5.1" (B4) | a new Section 2.5.2 |

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

*Script:* `analysis/calibration_fit_statistics.R` *Outputs:*
`output/calibration_fit_statistics.csv`,
`output/calibration_predictions.csv`,
`figures/calibration_observed_vs_predicted.png`, and
`app_data/calib_pred_by_variant.rds` (no longer 0 bytes)

The calibration-prediction step is restored and the statistics are now
reproducible from saved artifacts, without re-running `precompute.R`.

**The reproduction is exact.** Re-running the optimiser on the
reconstructed objective returns the saved parameters to machine
precision — SAR 0.0026852 and `rear_surv` 0.5427946, difference 0.00e+00
on both. So the numbers below are what this model actually does.

| Statistic, TDM-weighted, 2014–2024 | SI currently says | Regenerated |
|------------------------------------|-------------------|-------------|
| R² (squared correlation)           | 0.72              | **0.13**    |
| R² (Nash–Sutcliffe)                | —                 | **−0.59**   |
| RMSE                               | 8,240             | **13,345**  |
| MAPE                               | 24%               | **51%**     |

**The published values cannot be recovered under any variant tested.**
Not with all 14 years instead of 11 (R² = 0.31, RMSE = 11,829). Not
per-variant. And not by reverting the carrying capacity to the pre-rerun
12,493, which gives R² = 0.064, RMSE = 11,812, MAPE = 53%. They belong
to a structurally different earlier version of the pipeline, not to a
different parameterisation of this one.

A negative Nash–Sutcliffe means the model predicts observed escapement
**worse than the mean of the observations would**. The year-by-year
residuals show why:

| Year | Observed | Predicted | Residual    |
|------|----------|-----------|-------------|
| 2016 | 14,473   | 34,148    | +19,675     |
| 2023 | 37,321   | 17,169    | −20,152     |
| 2024 | 45,541   | 16,942    | **−28,599** |

The model produces a nearly flat 17,000–25,000 trajectory and misses the
2023–2024 surge entirely.

**What to do — pick one:**

1.  **Replace the sentence with the regenerated statistics and reframe
    the claim.** Honest and defensible if paired with the framing we
    already adopted for the projections: this is a *relative* comparison
    across alternatives under a fixed climatology, not a predictive
    hindcast. Two parameters are fitted to reproduce general magnitude;
    interannual variation is driven by ocean conditions the model does
    not represent. **This is my recommendation.**
2.  **Delete the fit statistics and the "captured the temporal dynamics"
    claim (B6) entirely**, and describe the calibration as a
    magnitude-matching exercise. Safest, loses the least defensible
    content.
3.  Do not quote R² = 0.72. It is not reproducible and the code is
    public.

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

## E. Corrections already applied on 2026-07-28 — verify, do not redo

Listed so it is clear what changed and by whom. **Everything in this
section is done.** The automated edits are highlighted **bright green**
in the documents; yellow and cyan highlighting predates them.

| Where | Change | By |
|----|----|----|
| Manuscript §3.2 | Composite scores → 0.573 / 0.534 / 0.530 / 0.502 / 0.500, and the full ranking | author |
| Manuscript, Fig. 3 caption | Martin threshold 12.8 → 12.14 °C | author |
| Manuscript §3.1.1 | Oct 16 – Nov 15 cooling benefit stated as −0.35 to −0.91 °C | author |
| Manuscript §3.1.2 | Martin-only index range added; PB1/PB3/PB6 falling below no-bypass explained by the front-loading mechanism | author |
| Manuscript, Discussion | EVPI range 0.000–0.027 | author |
| Manuscript §3.1 | Spearman rho 0.95 → **0.96** (recomputed on the current run, 0.9624) | **automated, green** |
| Manuscript, Results | EVPI point estimate qualified as the upper bound of the 0.000–0.027 range, so it stops contradicting the Discussion | **automated, green** |
| SI §S1 | SAR sentence replaced with verified brood-year values; bootstrap-CI note added; window corrected (data end at brood year 2019) | author |
| SI §S2.7/S2.8 | `rear_surv` 0.543 calibrated vs 0.5419 optimizer start clarified | author |
| SI, parameter table | Renumbered S2-7 → **S2-9** (S2-7 was used twice) | **automated, green** |
| SI, parameter table | Header "Calibrated Value" → "Value (empirical input)"; cells annotated with the calibrated estimates, which the prose already gave | **automated, green** |
| SI §S2.8 | `[EDITORIAL — DECISION NEEDED]` note inserted on the fit statistics | **automated, green** |
| Response | Crossovers 17.25/16.59 °C, the 12.2–16.6 °C band, operational range 14.1–18.2 °C | author |
| Response | Martin-only "functionally extinct" argument for the degenerate corner | author |
| Response | "late-spawning" → "November-spawning cohorts" throughout | author |
| Response | `[INTERNAL]` blocking list updated to current status | **automated, green** |

**Two verification notes.** The mixed-run problem is closed: no value
from a superseded run (12,493, 0.004161, 17,245, 18,396, 18,953, 0.537,
0.499, 0.491) survives anywhere in the manuscript or SI. And the single
remaining "late-spawning" string in the response sits inside an
`[INTERNAL, resolved]` note documenting the change — it is not a stale
error, and it goes when the `[INTERNAL]` notes are deleted (section C).

## F. Available and verified, but not currently used

Not errors. Each is a supported argument sitting unused, offered in case
a reviewer presses.

| \# | Argument | Source |
|----|----|----|
| F1 | **The projected-vs-observed gap is not uniform across TDM models.** Under Bratovich (0.51 of the panel weight) no-bypass projects 13,178, which is **59%** of the observed 2014–2024 mean of 22,491. The headline 34% is the model-averaged figure, dragged down by Martin projecting extinction under *every* alternative including no-bypass. Much easier to defend than the pooled number, and it pre-empts the obvious follow-up. | `output/frontloading_index_by_variant.csv` |
| F2 | **§3.1.1 could be strengthened.** It currently says late-November warming occurred "e.g., alternative PB1, PB3, and PB6", which is true as an example. In fact **all eight** bypass alternatives end warmer than no-bypass in late November, and the crossover date tracks each schedule's end date exactly — which is what makes the front-loading mechanism testable rather than anecdotal. | `output/frontloading_crossover_dates.csv` |
| F3 | **The cohort decomposition itself.** November cohorts hold 60.5% of redds but account for 97.5% of PB6's deficit relative to PB4; October cohorts are slightly *better* off. Currently asserted in the Discussion without the supporting numbers. | `output/frontloading_cohort_decomposition.csv`, `figures/frontloading_cohort_decomposition.png` |
| F4 | **EVPI normalisation caveat.** The choice of pooled vs per-state normalisation moves EVPI as much as the weight set does. Stating it pre-empts a fair criticism. | `output/evpi.csv` |

## G. Code and data items needing a decision before the repo goes public

The manuscript points at this repository, so these are now publication
decisions, not housekeeping.

| \# | Item | Status |
|----|----|----|
| G1 | **Alternative-specific spawn timing is computed and then discarded.** `SalmonCountR/precompute.R:659` splits redds by alternative and year; line 773 overwrites that with a split by year only, pooling all 36 alternatives. The CLM's temperature-driven shift in spawn timing therefore never affects any result. Flagged with a comment in the code; **not changed**, because correcting it would move every downstream number. Decide whether to disclose, fix, or defend it as holding spawn timing constant so differences are attributable to incubation exposure alone. | flagged in code |
| G2 | **Forecast temperatures are misaligned in leap years.** The forecast series is built by day-of-year, so leap years shift by one day against non-leap years — up to 1.08 °C on a given calendar date. Immaterial to the conclusions, but it is a real defect someone reading the code will find. | not fixed |
| G3 | **`sar_percent` in `app_data/SAR LAR Releases.xlsx` is identical to `sar`** — never multiplied by 100. Anyone reading that column as a percentage gets values 100× too small. | not fixed |
| G4 | **`app_data/*.rds` carry no provenance.** Nothing in the files records which data vintage or which commit produced them, which is the root cause of the mixed-run problem. See `APP_DATA_UPDATE_OPTIONS.md` §Option 2. | **plumbed 2026-07-29.** `refresh_data_year.R --apply` writes `app_data/data_vintage.rds`; `global.R` reads it and the app shows it in the footer, on the About tab, and at the head of every CSV export. Needs one applied refresh to populate — until then it reads "not recorded" |

*Fixed already:* `app_data/sim_redds.rds` was a stale artifact from an
older pipeline that did not reproduce `egg_summary.rds`; `precompute.R`
now writes `sim_redds` and `sim_future`. `calib_pred_by_variant.rds` is
no longer 0 bytes.

------------------------------------------------------------------------

## Where everything lives

|   | Path |
|----|----|
| Figures for the manuscript and SI | `figures/` |
| Tables and supporting numbers | `output/` |
| Scripts that regenerate them | `analysis/` |
| Full findings, organised as edits | `MANUSCRIPT_REVISION_HANDOFF.md` |
| Full findings, as an analysis record | `REVISION_FINDINGS.md` |
| Options memo for the app data question | `APP_DATA_UPDATE_OPTIONS.md` |
| Pre-edit copies of the three documents | `manuscript_backup/` |

**Note:** the live `.docx` files are no longer in the repo root — only
the pre-edit backups are. `.gitignore` excludes `*.docx`, so the working
copies live alongside the repo by design.

## Suggested order

1.  **Resolve D1 first.** It determines whether Figure S2-2 exists and
    what Section S2.8 says, and it is the only item that changes the
    scientific claims rather than the presentation.
2.  **A4, then A2.** Both are tables and both are copy-paste from files
    that already exist. Resolve D2 and D5 before Table S2-8 goes in. A4
    is the direct answer to the elicitation-rigour complaint. Fix the
    S2-9 cross-reference (B2) at the same time.
3.  **Table 3 and the six-callout renumbering** (current Table 3 becomes
    Table 4).
4.  **B1, A3, A1.** Swap Figures 3 and 5, insert Figure 6, and add the
    two callout paragraphs.
5.  **B6, B3, B4, C.** Cleanup, with B6 falling out of whatever D1
    decides, and the response placeholders filled last.

Once A1 through A4 are in, every `[N]` and `[S_]` in the response can be
filled with a real number, and the response stops making claims the
manuscript does not support.
