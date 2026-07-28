---
editor_options:
  markdown:
    wrap: 72
---

# Exhibit placement plan, 2026-07-28

Where each new figure and table goes, what number it gets, what its
caption says, and every callout that has to change. Numbering follows
order of first mention, which is what the journal copyeditor will check.

**Summary of numbering**

| Exhibit | Number | Effect on existing numbering |
|----|----|----|
| TDM model comparison table (A2) | **Table 3** (new) | Current Table 3 becomes **Table 4** |
| TDM weight sensitivity figure (A1) | **Figure 6** (new) | None |
| Replotted TDM curves (B1) | Figure 3, unchanged | None, but caption becomes two-panel |
| MCDA composite scores (A3) | Figure 5, unchanged | None |
| Panelist weights (A4) | Table S2-7, already numbered | None |
| Panelist justifications (A4) | Table S2-8, already numbered | None |
| Calibration observed vs predicted (D1) | **Figure S2-2** (new) | None |
| Front-loading decomposition (optional) | **Figure S2-1** (new) | None |

------------------------------------------------------------------------

## 1. Table 3 (new) — TDM model comparison

**Why Table 3 and not Table 4.** Tables are numbered by first mention.
Current order is Table 1 (§2.3), Table 2 (§2.4), Table 3 (§3.1). The
comparison table is first mentioned in §2.5.2, which precedes §3.1, so
it takes the number 3 and the consequence table becomes Table 4.

**Insert:** at the end of new Section 2.5.2 (Temperature-Dependent
Mortality Model Weighting), immediately after the paragraph describing
the Delphi allocation and the final weights. The table then sits
directly above the weights it explains, which is what Reviewer 2 asked
for.

**Source:** the table already drafted inside the response letter. Verify
the survival numbers against `output/figure3_etf_survival_by_temp.csv`.

**Add this callout sentence** at the end of that paragraph:

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

**Second callout, in the Discussion**, where the extrapolation argument
is made: append "(Table 3)" to the sentence identifying the 12.14 °C
threshold as a winter-run field estimate.

### Renumbering checklist — current Table 3 becomes Table 4

Four callouts in the manuscript and two in the SI:

| Document | Location | Current text |
|----|----|----|
| MS | §3.1, first paragraph | "…biological benefits and economic costs (Table 3)." |
| MS | §3.1.3 Hydropower Costs | "…scaled directly with bypass volume (Table 2, Table 3)…" |
| MS | §3.1.2 | "…worst to best alternative for Chinook salmon (Table 3)." |
| MS | §3.2 | "…little contribution from fish benefits (Table 3, Figure 5)." |
| SI | S2 | "…four meteorological scenarios as seen in Table 3 in the main text." |
| SI | S2 | "…in contrast to 7,600-11,073 seen in Table 3)." |

Plus the caption itself in the Tables section. If these were inserted as
Word cross-reference fields they will update on F9; if they are plain
text, change all six by hand.

------------------------------------------------------------------------

## 2. Figure 6 (new) — TDM weight sensitivity

**Why Figure 6.** First mention falls in §3.2, after Figure 5. No
renumbering.

**Insert:** in Section 3.2 (Trade-off Evaluation), immediately after the
paragraph reporting the composite scores and before the paragraph on
Reclamation's selection of PB4.

**File:** `figures/tdm_weight_sensitivity.png`. Supporting numbers in
`output/tdm_weight_sensitivity.csv` and
`output/tdm_weight_sensitivity_martin_sweep.csv`.

**Add this paragraph** to §3.2 as the callout:

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

**Third callout, in the Discussion**, alongside the existing robustness
sentence: "(Figure 6)".

------------------------------------------------------------------------

## 3. Figure 3 — replace image, rewrite caption

**No number change.** Swap the image for
`figures/figure3_tdm_curves.png` and replace the caption, which
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

**Existing callout in §2.5.3 needs no change**, but add a pointer to
panel (b) in the sentence describing the three formulations, since panel
(b) is the one that carries the argument.

------------------------------------------------------------------------

## 4. Figure 5 — replace image only

**No number or caption change.** Swap in
`figures/mcda_composite_scores.png`, which already carries the numeric
bar labels the caption promises and already reflects the current model
run. This also fixes the stale bar heights.

------------------------------------------------------------------------

## 5. SI Tables S2-7 and S2-8

Both are already referenced in the SI text; only the tables themselves
are missing.

**Insert both** in Section S2.5, immediately after Table S2-6 and after
the paragraph beginning "Individual weightings are reported in Table
S2-7."

**Sources:** `output/table_S2-7_panelist_weights.csv` and
`output/table_S2-8_justifications.md`.

**Captions:**

> *Table S2-7. Individual expert weightings of the three
> temperature-dependent mortality formulations. Weights were allocated
> independently in a single round, each summing to 1.0. The panel mean
> is the weighting used in the analysis.*

> *Table S2-8. Written justifications recorded by each panelist
> alongside their weight allocation, reproduced verbatim.*

**Before inserting S2-8, resolve D2 and D5.** Panelist 5's prose states
45/20/35 against a recorded row of 40/30/30, which becomes visible the
moment the justification appears next to the table. And decide whether
panelists are named or anonymized; if named, get consent, because the
justifications are candid.

------------------------------------------------------------------------

## 6. SI cross-reference fix (B2)

Section S2.8 currently reads "Final calibrated parameter values are
presented in Table S2-7." That number now belongs to the panelist
weighting table. Change to **Table S2-9**, then delete the `[EDITORIAL]`
note attached to the S2-9 caption.

------------------------------------------------------------------------

## 7. New SI figures

### Figure S2-1 (optional) — front-loading decomposition

Worth adding because the Discussion now asserts the front-loading
mechanism as a result rather than a conjecture, and nothing in the paper
shows it.

**Insert:** Section S2.4 or a short new subsection following the TDM
formulations. **File:** `figures/frontloading_cohort_decomposition.png`.

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

Panel (c) is the one that carries the argument; it shows the sign
reversal by spawn date and that December onward is exactly zero.

### Figure S2-2 — calibration observed versus predicted

Only if D1 option 1 is adopted, that is, reporting the regenerated
statistics rather than deleting them.

**Insert:** Section S2.8, with the revised fit statistics. **File:**
`figures/calibration_observed_vs_predicted.png`.

**Caption:**

> *Figure S2-2. Observed versus predicted in-river escapement over the
> 2014–2024 calibration period, TDM-weighted. The model reproduces the
> general magnitude of escapement but not its interannual variation;
> residuals are largest in 2016, 2023 and 2024.*

If D1 option 2 is chosen instead, skip this figure entirely, since a
plot showing a flat trajectory against a 2023–2024 surge invites exactly
the question the deletion is meant to avoid.

------------------------------------------------------------------------

## 8. Filling the placeholders in the response

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

------------------------------------------------------------------------

## Order of operations

1.  Resolve **D1** first. It determines whether Figure S2-2 exists and
    what Section S2.8 says, and it is the only item that changes the
    scientific claims rather than the presentation.
2.  Insert **Tables S2-7 and S2-8** (resolve D2 and D5 first), then fix
    the S2-9 cross-reference.
3.  Insert **Table 3** and run the six-callout renumbering.
4.  Swap **Figures 3 and 5**, insert **Figure 6**, add the two callout
    paragraphs.
5.  Fill the response placeholders and delete the `[INTERNAL]` and
    `[EDITORIAL]` notes.
