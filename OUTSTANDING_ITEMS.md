---
editor_options: 
  markdown: 
    wrap: 72
---

# What is still missing, updated 2026-07-29

Audited the three current files against every promise made in the
response to reviewers. The prose is in good shape. **The gap is almost
entirely figures and tables: we promise four new or revised exhibits to
the reviewers and none of them are in the documents yet.**

Confirmed by file inspection: the five images embedded in the manuscript
are byte-identical to the 2026-07-21 original, so no figure has been
replaced or added.

**Closed since the 07-28 pass and removed from this file:** D1 (fit
statistics — decided by deletion, so Figure S2-2 is no longer needed),
D3 (SAR provenance — reproducible, see B7 for the one error it exposed),
G4 (provenance stamp — built, needs one refresh run), and the itemised
07-28 correction list in section E. Anything still listed here is still
open.

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

**⚠️ SUPERSEDED 2026-08-18 — this paragraph changed qualitatively, not just
numerically.** Under the correction the sweep gains an intermediate regime that
did not exist before: PB1 holds the top rank to a Martin weight of **0.972**,
then **PB2c takes the top rank from 0.974 to 0.996**, and no-bypass only from
**0.998**. The Martin-only span also roughly doubles, to **3.5–171.1** spawners.
Replacement text:

> …Holding the other two formulations at their elicited ratio and increasing
> the weight on Martin et al. (2017) from 0 to 1, PB1 retained the top rank up
> to a weight of 0.972. Above that the ranking became unstable: PB2c held the
> top rank between 0.974 and 0.996, and the no-bypass alternative only from
> 0.998. Across that upper range the projected adult population index spans
> just 3.5 to 171 spawners across all nine alternatives, so the Chinook
> objective is normalized across a range of a few hundred fish against a
> baseline of tens of thousands, and the ranking is decided almost entirely by
> the hydropower objective.

Source: `output/tdm_weight_sensitivity_martin_sweep.csv`. Note the conclusion —
the elicited weighting is nowhere near the unstable region — is unchanged and if
anything better supported, since the instability now shows up as two regime
changes rather than one.

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
| B7 | **The SAR paragraph in SI §S2.6 mixes two brood-year windows.** The group count, mean and standard deviation are over brood years 2008–2019 (n = 27, mean 0.25%, sd 0.0024) — all three reproduce exactly. The stated *range*, "0.02% (brood year 2017) to 0.56% (brood year 2016)", is the range over brood years **2011–2019 only**. Over the window the sentence actually claims, the maximum is **0.95% in brood year 2010**; BY2016 is second. **Fix:** change the range to 0.02% (BY2017) to 0.95% (BY2010) and keep the window. Restricting the window instead would make n = 24 and the mean 0.22%, which breaks the sentence tying the mean to the 0.0025 optimizer starting value. | `analysis/sar_from_cwt.R`, `output/sar_by_brood.csv` |

## C. Placeholders still in the text

| Where | What |
|----|----|
| Response | `[N]` ×3 and `[S_]` ×2 — resolve once the figures and tables are in |
| Response | `[Consequence Models]` — remove the brackets, the rename is done |
| Response | All `[INTERNAL: ...]` notes — delete before submission |
| SI | `[citation to be added]` for the Hazel Avenue statistics → the WTMP deliverable |
| SI | `[EDITORIAL: ...]` note on Table S2-9 numbering — delete once B2 is resolved |
| SI | `[EDITORIAL — DECISION NEEDED]` note on the fit statistics — delete; D1 is decided |
| Manuscript | `[*add repository URL upon acceptance*]` — normal, leave until acceptance |

## D. Decisions that need a human

### D1. Calibration fit statistics — **DECIDED 2026-07-29: deleted**

The author reports having deleted the fit statistics from the SI, and
describing the calibration as magnitude-matching instead. Defensible:
the published values were not reproducible under any variant tested,
which is consistent with their coming from a structurally earlier
pipeline.

**All that is left is confirmation.** Only the pre-edit backups are in
the repo, so none of the following could be verified here — three things
to check in the live SI:

1.  The fit-statistics sentence is gone (R² 0.72, RMSE 8,240, MAPE 24%).
2.  **B6 is gone with it** — the separate claim that the model "captured
    the magnitude and general temporal dynamics… including major peaks
    in 2015, 2018, and 2023–2024". It is a different sentence and
    deleting the statistics does not remove it.
3.  The `[EDITORIAL — DECISION NEEDED]` note on §S2.8 is deleted.

**Keep these numbers to hand** in case a reviewer asks why no fit
statistics are reported. Regenerated from
`analysis/calibration_fit_statistics.R`, TDM-weighted, 2014–2024: R²
0.13 (published 0.72), Nash–Sutcliffe **−0.59**, RMSE 13,345 (published
8,240), MAPE 51% (published 24%). A negative Nash–Sutcliffe means the
model predicts observed escapement worse than the mean of the
observations would; it produces a nearly flat 17,000–25,000 trajectory
and misses the 2023–2024 surge entirely (2024 under by 28,599). The
published values are not recoverable under any variant tested, so the
defensible framing is magnitude-matching, not hindcast skill.

### D2, D4–D7

| \# | Item |
|----|----|
| D2 | **Panelist 5's weights**: prose says 45/20/35, recorded row says 40/30/30. Visible the moment Table S2-8 reproduces the justification verbatim. Source: `output/table_S2-8_justifications.md` |
| D4 | **Acknowledgments** omit one panel member and spell another differently from the scoresheet. |
| D5 | **Anonymise or name the panelists** in Tables S2-7/S2-8. If named, we need consent; the justifications are candid. |
| D6 | **Three highlight colours** now layer in the manuscript (yellow, cyan, green) with no key. Green marks the 2026-07-28 automated edits. Collapse or add a legend before this goes to Brian or the editor. |
| D7 | **Title scoping.** Reviewer 1 objects that "Balancing" implies absolute performance. Conclusions were scoped; the title was left alone deliberately. |

------------------------------------------------------------------------

## E. Already applied on 2026-07-28 — do not redo

Sixteen numeric and wording corrections went into the three documents on
2026-07-28: composite scores and ranking, the 12.14 °C threshold, the
cooling-benefit and EVPI ranges, Spearman rho 0.96, the SAR sentence,
the `rear_surv` clarification, the S2-7 → **S2-9** renumbering and
parameter table header, and the crossover temperatures and
"November-spawning" wording in the response. The itemised list is in
this file's git history at commit `08a7f4b`; the automated edits are
highlighted **bright green** in the documents, and yellow and cyan
predate them (D6).

**Two verification notes still worth having.** The mixed-run problem is
closed: no value from a superseded run (12,493, 0.004161, 17,245,
18,396, 18,953, 0.537, 0.499, 0.491) survives anywhere in the manuscript
or SI. And the single remaining "late-spawning" string in the response
sits inside an `[INTERNAL, resolved]` note documenting the change — it
is not a stale error, and it goes when the `[INTERNAL]` notes are
deleted (section C).

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
| G1 | ~~**Alternative-specific spawn timing is computed and then discarded.**~~ **CLOSED 2026-08-18 — corrected.** Redds were pooled across all 36 alternatives before survival was computed, so the CLM's temperature-driven shift in spawn timing never reached any result. Alternative-specific spawn timing is now the default; `ARG_G1_ALT_SPAWN=0` reproduces the superseded numbers. Measured across five seeds in both arms: PB3−PB5 reverses sign (+43 → −261), PB6/PB4 becomes 9,199/9,550, PB2b−PB5 reverses in the composite, efficiency becomes 45/79/84, Spearman ρ 0.9624 → 0.9205. PB4 (the selection) moves +3.7 adults; PB1 still ranks first. Replacement text in `G1_DISCLOSURE.md`, numbers in `G1_FINDINGS.md` and `output/g1_revision_*.csv`. **Document changes this forces are listed in §H.** | `G1_FINDINGS.md`, `analysis/g1_revision_numbers.R` |
| G1a | ~~**The cohort decomposition lost its exactness.**~~ **CLOSED 2026-08-18 — recomputed exactly.** The PB6−PB4 split was exact only because alternatives shared one redd distribution; G1 removed that. `frontloading_cohort_decomposition.R` Part 3b now computes the exact two-channel identity (asserted in code) and it changed the finding materially: composition is **41%** of the PB6−PB4 gap, and the December cohorts went from 0% to **47%** of it — their survival is untouched but front-loading moves redds out of the safe December window. The old 49%/49% November figures are superseded. Numbers in `output/frontloading_cohort_two_channel.csv`; replacement table and mechanism sentence in `MANUSCRIPT_REVISION_HANDOFF.md` §2.4. | `analysis/frontloading_cohort_decomposition.R` Part 3b |
| G2 | **Forecast temperatures are misaligned in leap years.** The forecast series is built by day-of-year, so leap years shift by one day against non-leap years — up to 1.08 °C on a given calendar date. Immaterial to the conclusions, but it is a real defect someone reading the code will find. | not fixed |
| G3 | **`sar_percent` in `app_data/SAR LAR Releases.xlsx` is identical to `sar`** — never multiplied by 100. Anyone reading that column as a percentage gets values 100× too small. **No live consumer is affected:** `analysis/sar_from_cwt.R:32` drops the column outright (`select(-any_of("sar_percent"))`) with the reason in a comment, recomputes the percentage itself, and `analysis/data_sources.R:98` records the defect in the source register. The remaining exposure is a human opening the workbook and trusting the header. **Fix when the workbook is next regenerated:** either multiply by 100 or delete the column — do not edit the snapshot in place, since it is a provenance-stamped input. | latent, defended in code |
| G4 | ~~**`app_data/*.rds` provenance**~~ — **DONE 2026-08-14.** `refresh_data_year.R --apply` was run and wrote `app_data/data_vintage.rds` (8 sources, stamped at commit `7b064a2`); the app now shows a real vintage in the footer, on About, and atop every CSV export instead of "not recorded". The run replaced no snapshots — nothing new was available to download — so no number moved. Report: `output/data_refresh_report_2026-08-14.md`. | **closed** |

------------------------------------------------------------------------

## H. Document changes forced by the G1 correction — **all still open**

The code side of G1 is finished and verified. Everything in this section is
Word-document work that nobody has done yet. Numbers are means over five seeds;
regenerate any of them with `analysis/g1_revision_numbers.R`.

Two sanity anchors before editing: running `ARG_G1_ALT_SPAWN=0` reproduces the
**submitted** values exactly (Table 3 for all nine alternatives, Spearman
ρ = 0.9624), and `main` is tagged `as-submitted-2026-07-28`. So every "was"
below is verifiably the number now in the documents.

### H0. Reporting convention — **DECIDED 2026-08-18**

**Levels come from the single committed run. Spread is quoted alongside.
Differences are reported as paired contrasts, not as subtractions of levels.**

Generate every level with `analysis/reporting_values.R`, which reads the
committed `app_data` and writes `output/reporting_values.csv`. It runs in
seconds and needs no seed snapshots, so anyone who clones the repo reproduces
the paper's tables directly.

| Quantity | Source | Why |
|---|---|---|
| Table 3, composites, efficiency, ranges | committed run, via `reporting_values.R` | table, every figure, the app and the repo all trace to one artefact |
| Uncertainty on a level | ±173–217 fish, from the five-seed replication | stating the range is what keeps a single run honest |
| PB3−PB5, PB6−PB4, and any interpreted difference | **paired within-seed mean** from `g1_revision_numbers.R` | see below |

**Why differences are handled separately.** The seed shift is largely *common*
across alternatives, so it cancels in a paired difference. Levels carry ±173–217
fish of run-to-run noise; the contrasts that carry the argument carry only
±21–36. Subtracting two numbers out of Table 3 throws that cancellation away and
would understate what the model can actually resolve. Report PB5 leading PB3 by
**261 ± 36**, not by the 240 you get by subtracting the table.

**Why not the five-seed mean for levels.** It was the other serious candidate.
Against it: the paper's numbers would then correspond to no single committed
artefact, so reproducing Table 3 would mean a ~70-minute five-seed rebuild rather
than opening `app_data`; and Figure 5 and the composite scores would need
seed-averaging machinery that does not exist (`mcda.R` reads `app_data`
directly). Both headline claims hold at the committed run anyway — PB3−PB5 is
−240 there against a five-seed mean of −261, PB6−PB4 is −335 against −352, both
the right sign and inside the interval — so averaging buys robustness the
argument does not need.

**Put the five-seed replication table in the SI** as the uncertainty
characterisation. That answers a reviewer objecting to a single realisation at no
cost to internal consistency.

*Not a factor in this decision:* the Shiny app diverging from the paper. The app
continues to be developed past the journal version, so the two are expected to
drift.

| \# | Location | Was | Becomes | Kind of change |
|----|----|----|----|----|
| H1 | MS Table 3, abundance column | 7,600 / 8,560 / 10,505 / 10,974 / 11,073 / 9,116 / 9,545 / 9,080 / 9,350 (NB→PB6) | **7,412 / 8,306 / 10,443 / 10,992 / 11,274 / 8,954 / 9,427 / 9,194 / 9,091** | retype all nine; add the ±173–217 run-to-run note |
| H1a | MS §3.2, composite scores | PB1 0.573, PB2 0.534, PB4 0.530 (values from the last pass) | PB1 **0.555**, PB4 **0.514**, PB2 **0.513**, PB3 0.508, PB2c 0.502, NB 0.500, PB5 0.479, PB2b 0.465, PB6 0.403 | retype |
| H1b | MS §3.2, MCDA ordering | PB1 > PB2 > PB4 > PB3 > PB2c > NB > PB2b > PB5 > PB6 | PB1 > **PB4 > PB2** > PB3 > PB2c > NB > **PB5 > PB2b** > PB6 | **two swaps** — see caveat below |
| H2 | MS Discussion, PB3 vs PB5 | "differed by fewer than 40 adults… consequential at specific points" | PB5 leads by 261 ± 36; schedule matters throughout | **rewrite, sign reversed** |
| H3 | MS Discussion, PB6 vs PB4 | "9,350 against 9,545" | "**9,091 against 9,427**"; the gap itself is **−352 ± 21** (paired, per H0) | retype |
| H4 | MS Discussion, efficiency | "roughly 47… against 76 for PB4 and PB2c" | PB6 **45**, PB4 **78**, PB2c **84** | **restructure** — PB4 and PB2c no longer share a value |
| H5 | MS Discussion, volume correlation | "Spearman ρ = 0.95" (edited to 0.96 last pass) | 0.92 | retype |
| H6 | MS Discussion, mechanism | one channel (worse incubation window) | two channels: worse window **and** redds displaced into it, ~59% / 41% | **rewrite** — see handoff §2.4 |
| H7 | MS + SI, "late-spawning cohorts" | attributed to late spawners | November window, but say *which channel* | rewrite; the old justification that December spawners are unaffected is **false** |
| H8 | SI, cohort decomposition table | Nov 1–15 49%, Nov 16–30 49%, Dec 0% | Nov 16–30 44%, Dec 1–15 32%, Dec 16–Jan 15%, Nov 1–15 8% | replace table |
| H9 | MS Methods, spawn timing | pooling not described | state that each alternative uses its own simulated redds | **new sentence** — draft in `G1_DISCLOSURE.md` §1 |
| H10 | Response to reviewers | — | disclose the correction | **new section** — draft in `G1_DISCLOSURE.md` §2 |
| H11 | SI reproducibility | — | `ARG_G1_ALT_SPAWN` / `ARG_SEED`, five-seed reporting | **new note** — draft in `G1_DISCLOSURE.md` §3 |
| H12 | Anywhere quoting a difference between alternatives | subtraction of two table values | **paired within-seed contrast with its own spread** (per H0) | **differences under ~400 fish are not supportable from a subtraction of levels** — true of the submitted numbers too |
| H13 | All figures embedded in MS/SI | 2026-07-21 images | regenerated versions in `figures/` | re-embed; see §A, none of the promised exhibits are in the documents yet either |
| H14 | MS §4.3 + Response, EVPI | range 0.000–0.027, upper bound 0.026 (4.6%) | range **0.000–0.034**, upper bound **0.034 (6.1%)**; two of four combinations now give zero, not one | retype + adjust the sentence about how many cases are zero |
| H15 | Figure 5 / MCDA composite chart | pre-correction bar values | regenerated `figures/mcda_composite_scores.png` | already regenerated in the repo; re-embed |
| H16 | MS §3.2, Figure 6 sensitivity paragraph | PB1 to 0.988, then NB at 1.000; span 3.6–84.8 | PB1 to **0.972**, **PB2c 0.974–0.996**, NB from **0.998**; span **3.5–171** | **rewrite** — a regime that did not exist appears; replacement text in §A-bis above |
| H17 | MS §2.3 / SI, Martin index vs net hazard | "exact inverse order, ρ = −1.00" | ρ = **−0.983**; PB1 and PB3 swap | **soften the claim** — it is now falsifiable against the repo as written |
| H18 | SI S2, "in contrast to 7,600–11,073 seen in Table 3" | 7,600–11,073 | **7,412–11,274** | retype |
| H19 | Anywhere quoting projected vs observed | NB Bratovich 13,178 = 59% of observed; model-averaged = 34% | **12,828 = 57%**; **33%** | retype |

**H12 is the one with teeth.** Pooling suppressed run-to-run noise about
fourteen-fold (±15–21 fish reported, ±173–217 actual). Any comparison in the
current text resting on a margin under ~400 fish was never resolvable from a
single run, and reporting it to the digit overstates precision. This applies to
the submitted numbers as much as the corrected ones.

**Caveat on H1b — only one of the two swaps is real.** PB5 overtaking PB2b is a
genuine, reproducible effect of the correction: same direction at all five seeds,
3.5–10.9× the run-to-run noise. PB4 overtaking PB2 is **not** — they differ by
0.0009 at the committed run and 0.0002 in the five-seed mean, against a composite
standard deviation of 0.002–0.003. Their order is decided by the seed, not by the
alternatives. Do not describe PB4 as ranking above PB2. Either present them as
tied, or report the ordering only to the resolution the model supports.

Note also that NB and PB2c are pinned at 0.500 and 0.502 **by construction** —
they are the min–max normalisation anchors, so the composite understates their
sensitivity and their apparent stability is not evidence of anything.

------------------------------------------------------------------------

## Where everything lives

|                                        | Path                             |
|----------------------------------------|----------------------------------|
| Figures for the manuscript and SI      | `figures/`                       |
| Tables and supporting numbers          | `output/`                        |
| Scripts that regenerate them           | `analysis/`                      |
| Full findings, organised as edits      | `MANUSCRIPT_REVISION_HANDOFF.md` |
| Full findings, as an analysis record   | `REVISION_FINDINGS.md`           |
| Options memo for the app data question | `APP_DATA_UPDATE_OPTIONS.md`     |
| Pre-edit copies of the three documents | `manuscript_backup/`             |

**Note:** the live `.docx` files are no longer in the repo root — only
the pre-edit backups are. `.gitignore` excludes `*.docx`, so the working
copies live alongside the repo by design.

## Suggested order

1.  **A4, then A2.** Both are tables and both are copy-paste from files
    that already exist. Resolve D2 and D5 before Table S2-8 goes in. A4
    is the direct answer to the elicitation-rigour complaint. Fix the
    S2-9 cross-reference (B2) at the same time.
2.  **Table 3 and the six-callout renumbering** (current Table 3 becomes
    Table 4).
3.  **B1, A3, A1.** Swap Figures 3 and 5, insert Figure 6, and add the
    two callout paragraphs.
4.  **B6 and B7, then B3, B4, C.** B6 is the D1 confirmation; B7 is a
    one-number fix in the SAR paragraph. Response placeholders last.

Once A1 through A4 are in, every `[N]` and `[S_]` in the response can be
filled with a real number, and the response stops making claims the
manuscript does not support.
