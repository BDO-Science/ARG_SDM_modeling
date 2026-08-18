# Folsom Bypass revision — analysis findings

**Manuscript:** "Balancing salmon conservation and hydropower", major revision at
*River Research and Applications*
**Prepared:** 2026-07-28
**Source of truth:** `ARG_SDM_modeling` repo at commit `26d61b8` ("Update PHABSIM K
source, rerun precompute, and improve figure readability")

This document reports what the model output actually says, task by task, against
the claims currently in the manuscript, SI and response to reviewers. Every
number below was recomputed from `SalmonCountR/app_data/*.rds` unless marked
otherwise. Scripts that regenerate each result are named per section.

> ## ⚠️ Read this before quoting any number below
>
> **This is a record of the 2026-07-28 pass and it predates the G1 correction
> (2026-08-18).** Commit `26d61b8` — the "source of truth" named above — pooled
> the simulated redd set across all 36 alternatives before computing egg
> survival, which severed the temperature → spawn-timing → mortality channel.
> That was a defect, and it has since been fixed.
>
> **Every model-derived number in this document is therefore the pre-correction
> value.** The findings and the reasoning still stand — that is why this file is
> kept — but the figures attached to them have moved, several materially and two
> qualitatively (the PB3/PB5 comparison reverses sign; the cohort decomposition
> gained a second channel).
>
> **`OUTSTANDING_ITEMS.md` §H is the authoritative list of current values.**
> Where it and this document disagree, §H wins. `G1_FINDINGS.md` explains what
> changed and why. To reproduce the numbers as written below, run with
> `ARG_G1_ALT_SPAWN=0`, or check out the tag `as-submitted-2026-07-28`.

**Status: all twelve tasks complete.** For the version organised as manuscript
edits rather than as an analysis record — which is the one to hand to someone
making the revisions — see `MANUSCRIPT_REVISION_HANDOFF.md`.

---

## ⚠️ READ THIS FIRST — the manuscript mixes two model runs

This is the single most important finding and it affects several sections at once.

The repo has been re-run twice since the MCDA numbers in Section 3.2 were
written:

| Commit | What changed | K (spawners) | Calibrated SAR | NB adult index |
|---|---|---|---|---|
| `c658df3` / `e0cd8c8` | *(the run Section 3.2 was written from)* | 12,493 | 0.004399 | 18,396 |
| `16e7a70` | ATU boundaries → SacPAS v3 (400 / 558 / 958) | 12,493 | 0.004161 | 17,245 |
| `26d61b8` | PHABSIM K source updated, precompute re-run | **33,185** | **0.002685** | **7,602** |

Consequence: **the consequence-table abundances (7,600–11,073) are from the
current run, but the Section 3.2 composite scores (0.537 / 0.502 / 0.500 / 0.499
/ 0.491) are from the pre-rerun run.** The two are not mutually consistent and a
reviewer comparing them will see it.

Everything in this document is computed from the **current** run (`26d61b8`).
Adopt these numbers throughout, or re-run and re-derive everything from a single
commit — but do not leave the two mixed.

---

## Task 1 — TDM weight sensitivity (highest priority)

*Script:* `analysis/tdm_weight_sensitivity.R`
*Outputs:* `figures/tdm_weight_sensitivity.png`, `output/tdm_weight_sensitivity.csv`,
`output/tdm_weight_sensitivity_martin_sweep.csv`

Objective weights held at 0.40 salmon / 0.50 hydropower / 0.10 steelhead;
meteorological years at 0.25 each. Composite score by alternative:

| Alt | Bratovich only (1,0,0) | Bartholow only (0,1,0) | **Martin only (0,0,1)** | Equal (⅓,⅓,⅓) | Elicited (.51,.24,.25) |
|---|---|---|---|---|---|
| NB | 0.5000 | 0.5000 | **0.5153** | 0.5000 | 0.5000 |
| PB1 | **0.6191** | **0.5122** | 0.4671 | **0.5525** | **0.5734** |
| PB2 | 0.5659 | 0.4825 | 0.2868 | 0.5157 | 0.5339 |
| PB2b | 0.4938 | 0.4537 | 0.3203 | 0.4722 | 0.4824 |
| PB2c | 0.4885 | 0.5017 | 0.5017 | 0.5017 | 0.5017 |
| PB3 | 0.5786 | 0.4468 | 0.3496 | 0.4967 | 0.5227 |
| PB4 | 0.5741 | 0.4664 | 0.3271 | 0.5077 | 0.5296 |
| PB5 | 0.4807 | 0.4389 | 0.3266 | 0.4555 | 0.4646 |
| PB6 | 0.5061 | 0.3291 | 0.2290 | 0.3960 | 0.4306 |

### The single number that matters

**Yes — the top-ranked alternative DOES change under 100% Martin.** It flips from
PB1 to **NB (0.5153)**, with PB2c second at 0.5017. Margin 0.0135.

So the composite ranking is *not* unconditionally invariant. But the honest and
still-strong version of the claim is:

> Holding Bratovich : Bartholow at the elicited 0.51 : 0.24 and sweeping the
> Martin weight from 0 to 1, **PB1 remains top-ranked for Martin weight 0 through
> 0.988.** PB2c takes the top rank over 0.990–0.998, and NB only at 1.000
> exactly. The elicited Martin weight is 0.25, where PB1 leads the runner-up by
> 0.0395.

**Suggested opening line for the response:** the preferred alternative is
invariant to the TDM weighting over 98.8% of the weight simplex along the
elicited axis, and changes only in the degenerate corner where the panel's other
two models receive zero weight — a weighting no panelist proposed (the lowest
individual Bratovich+Bartholow allocation was 0.60).

### ⚠️ Added after Task 12: what the 100%-Martin corner actually contains

*Source:* `analysis/frontloading_cohort_decomposition.R` Part 5 →
`output/frontloading_index_by_variant.csv`

The adult index, median of the final 20 forecast years, broken out by TDM model:

| Alt | Bratovich | Bartholow | **Martin** | Model-averaged |
|---|---|---|---|---|
| PB2c | 16,909 | 10,117 | **84.8** | 11,073 |
| PB2b | 17,037 | 9,470 | **49.6** | 10,974 |
| PB2 | 16,714 | 8,229 | **21.4** | 10,505 |
| PB4 | 15,769 | 6,254 | **8.0** | 9,545 |
| PB6 | 15,851 | 5,272 | **3.6** | 9,350 |
| PB3 | 15,401 | 5,250 | **4.0** | 9,116 |
| PB5 | 14,978 | 5,994 | **10.2** | 9,080 |
| PB1 | 14,686 | 4,454 | **4.5** | 8,560 |
| NB | 13,178 | 3,655 | **6.7** | 7,600 |

The model-averaged column reproduces the consequence table exactly (7,600 →
11,073), which independently validates those numbers.

**Under the Martin model alone every alternative is functionally extinct.** The
entire nine-alternative range is 3.6 to 84.8 spawners, and the trajectory never
equilibrates — NB declines monotonically across the whole 114-year run, reaching
about 7 spawners at the end of it.

This reframes the Task 1 result. At 100% Martin weight the Chinook objective is a
min-max normalisation across a range of a few dozen fish, so it carries almost no
information; NB scores 0.038 on that normalised axis and wins the composite purely
on the hydropower objective, where it scores 1.000 by construction. **The corner
where the preferred alternative flips is not a corner where no-bypass is better
for salmon — it is a corner where the salmon axis has collapsed and hydropower
decides.** That is a considerably stronger argument for the response than "no
panelist proposed that weighting," and both points can be made together.

---

## Task 2 — Section 3.2 composite scores

*Script:* `analysis/mcda.R` (bar labels added)
*Output:* `figures/mcda_composite_scores.png`, `output/mcda_composite_scores_v2.csv`

**The text is stale; the figure is fine.** This reverses the working assumption.

Quoted in Section 3.2 vs the run it came from vs current:

| Alt | Text | `e0cd8c8` run | **Current run** |
|---|---|---|---|
| PB1 | 0.537 | 0.541 | **0.5734** |
| PB2c | 0.502 | 0.502 ✓ | **0.5017** |
| NB | 0.500 | 0.500 ✓ | **0.5000** |
| PB2 | 0.499 | 0.499 ✓ | **0.5339** |
| PB4 | 0.491 | 0.491 ✓ | **0.5296** |

Four of the five match `e0cd8c8` to three decimals, so that is unambiguously the
source run; PB1's 0.537 vs 0.541 is a transcription slip.

**Full current ranking:** PB1 0.573 > PB2 0.534 > PB4 0.530 > PB3 0.523 >
PB2c 0.502 > NB 0.500 > PB2b 0.482 > PB5 0.465 > PB6 0.431.

Two things follow:
1. The old top-five really were separated by ~0.003 and genuinely could not be
   read off a bar chart — which is why the "figure can't resolve it" reading felt
   right. That is no longer the situation: **PB1 now leads PB2 by 0.040**, easily
   visible.
2. If Reviewer 1 read a different order off Figure 5, they may have been reading a
   correctly regenerated figure against stale text. Brian's "Needs to be fixed"
   is correct, and the fix is to the text.

Numeric labels are now drawn above each bar, as the revised caption promises.

---

## Task 3 — Hazel Avenue regression fit statistics

**Not in the repo. Cite the source instead.**

There is no monthly multiple regression for Hazel Avenue *or* Watt Avenue
anywhere in: the `ARG_SDM_modeling` repo, the sibling copy at
`Documents/R/reclamation/gitlab/american_river_SDM`, or any project spreadsheet
(`ARG_LAR_TempModeling_*.xlsx`, `SDM Power Bypass Temperature Modeling
Results.xlsx` — all contain scenario temperature series only, no regression
sheets).

The repo consumes NWIS gauge observations (station 11446500 Hazel, 11446980
Watt) and CE-QUAL-W2 scenario output; it never fits the regression. The Watt
Avenue R² = 0.94 / RMSE = 0.76 °C / MAE = 0.58 °C must therefore come from the
WTMP deliverable.

**Action:** replace the `[TBD]` placeholders in SI S1 with a citation to the WTMP
deliverable, or request the Hazel statistics from whoever produced it. Do not
compute a substitute — it would not be the same model.

---

## Task 4 — `rear_surv`: 0.543 or 0.5419

**Your edit was correct. 0.543 is the calibrated value; 0.5419 is the optimizer
starting value.**

Evidence, `SalmonCountR/precompute.R:983-991`:

```r
opt_combined <- optim(
  par     = c(0.0025, 0.5419),      # Starting values for SAR and rear_surv
  fn      = combined_sse, ...
  lower   = c(0.0001, 0.01),
  upper   = c(0.1, 1.0)
)
```

`calib_results.rds` holds the results: **SAR_mean = 0.002685214, rear_surv =
0.5427946** — the latter rounding to 0.543. Both values are used for all three
TDM variants.

Keep S2.7 as you edited it. Correct S2.8 if it disagrees.

---

## Task 5 — Table S2-7 SAR row

**Solved.** The source is `SalmonCountR/app_data/SAR LAR Releases.xlsx`.

Filtering to the **27 AMERICAN R release groups** gives mean = **0.002501** and
SD = **0.002371** — exactly the model's `SAR_mean` 0.0025 and `SAR_sd` 0.00237.

The **0.0018–0.0034 range is a 95% bootstrap percentile CI on the mean**
(reproduces as 0.00168–0.00341; mean ± 1.96 SE gives 0.00160–0.00339). It is an
uncertainty interval on a mean, *not* a spread of observations, so it is not
comparable to an interannual range and the table must say which it is.

**Recommended table structure (as you proposed — confirmed correct):**

| Parameter | Value | Basis |
|---|---|---|
| SAR, input | 0.0025 (SD 0.00237) | mean of 27 American River CWT release groups |
| SAR, 95% CI on mean | 0.0018–0.0034 | bootstrap percentile CI |
| SAR, calibrated | 0.00269 | optimizer result |
| `rear_surv`, input | 0.5419 | optimizer starting value |
| `rear_surv`, calibrated | 0.543 | optimizer result |

**⚠️ The 0.05%–0.68% range is not in this file, and two other things in the same
SI sentence are also wrong.** The sentence reads: *"The 2011–2024 mean SAR was
0.25% (0.0025), with substantial interannual variation (range: 0.05% to 0.68%)."*

Twelve aggregations were tried — release-group min/max, Q1–Q3, 5th–95th,
10th–90th, brood-year means, pooled brood-year SAR (returns ÷ releases), with and
without Discovery Park, over several brood windows, and mean ± 1 SD. None
produces 0.05%–0.68%. The closest are 0.02%–0.56% (pooled brood-year) and
0.02%–0.62% (release groups); the published range is shifted up at both ends. The
mean reproduces exactly, so the file is the right source and the range came from
somewhere else.

Separately, **the "2011–2024" window does not match the data.** The file covers
brood years 2008–2019 (releases 2009–2020), and the 0.0025 mean is over all 27
American River groups *including* broods 2008–2010. There is no 2020–2024 CWT
data in it.

**Verified replacement:** mean 0.25% (SD 0.0024) across 27 American River release
groups, brood years 2008–2019; pooled brood-year SARs range 0.02% (BY2017) to
0.56% (BY2016). Full brood-year table in
`MANUSCRIPT_REVISION_HANDOFF.md` §3.3.

**Spreadsheet bug:** `sar_percent` in `SAR LAR Releases.xlsx` is identical to
`sar` — never multiplied by 100.

---

## Task 6 — Elicitation appendix tables

*Script:* `analysis/elicitation_tables.R`
*Outputs:* `output/table_S2-7_panelist_weights.csv`, `output/table_S2-8_justifications.md`

### Table S2-7 (anonymized)

| Panelist | Martin (2017) | Bratovich (2020) | Bartholow & Heasley (2006) | Sum |
|---|---|---|---|---|
| Panelist 1 | 0.25 | 0.25 | 0.50 | 1.00 |
| Panelist 2 | 0.10 | 0.90 | 0.00 | 1.00 |
| Panelist 3 | 0.20 | 0.70 | 0.10 | 1.00 |
| Panelist 4 | 0.30 | 0.40 | 0.30 | 1.00 |
| Panelist 5 | 0.40 | 0.30 | 0.30 | 1.00 |
| **Panel mean (weights used)** | **0.25** | **0.51** | **0.24** | 1.00 |

Column means confirmed. All five allocations sum to 1.

### Two protocol problems

**1. SI Section S2.5 does not describe what was done.** There is no
criteria-scoring sheet and no Round 2 sheet anywhere in the manuscript project
folder or the repo — only `ScoreSheet_ARG_TDM_Round1.xlsx`, which contains direct
0–1 allocations. So S2.5's description of panelists scoring each model 0–100
against four criteria and then normalizing, followed by a revision round after
anonymous feedback, must be corrected to: **a single round of direct 0–1 weight
allocation with written justification.** Drop the revision step.

**2. Panelist 5's numbers and prose disagree.** Their justification states the
Martin model "received the highest weight (45%)", Bratovich "20%" and Bartholow
"35%", but the recorded row is 0.40 / 0.30 / 0.30. If the prose values were used
the panel mean would be **0.26 / 0.49 / 0.25** instead of 0.25 / 0.51 / 0.24.
The effect on results is small but the discrepancy is visible in Table S2-8
if the justification is reproduced verbatim. Decide which is authoritative and
either correct the sheet or add an editorial note.

---

## Task 7 — Egg-to-fry survival comparison

**Reproduces from the implementation.** Recomputed by calling
`SalmonCountR/functions.R::tdm_exp()` and `tdm_lin_martin()` directly with
ATU-paced stage boundaries (hatch 400, emergence 958), not from the SI parameter
list:

| T (°C) | Bratovich | Bartholow | Martin |
|---|---|---|---|
| 13 | 95.0 | 94.9 | 19.1 |
| 14 | 84.7 | 81.8 | 3.6 |
| 15 | 59.0 | 46.2 | 0.9 |
| 16 | 18.0 | 5.0 | 0.2 |
| 17 | 0.4 | 0.0 | 0.1 |

Maximum discrepancy against the table in the response: **0.24 percentage points**.
The table can stand as written.

**Crossovers confirmed:** Bratovich falls below Martin at **17.25 °C**, Bartholow
at **16.59 °C** (response says ~17.2 and ~16.6). This is load-bearing against
Reviewer 2's "essentially no temperature response" claim and it holds.

One refinement for precision: there is also a *lower* crossover at 12.15 °C,
below which Martin predicts 100% survival while the exponential models are
already below 100%. The correct phrasing is that **Martin is harsher than the
exponential models between roughly 12.2 °C and 16.6 / 17.3 °C, and gentler
outside that band.**

---

## Task 8 — Figure 3 replotted

*Script:* `analysis/figure3_tdm_curves.R`
*Outputs:* `figures/figure3_tdm_curves.png`, `output/figure3_etf_survival_by_temp.csv`

- **The plotting code already used 12.14**, matching
  `functions.R::tdm_lin_martin`. The caption's 12.8 was the error; your
  correction to 12.14 is right.
- x range extended to **10–18 °C**.
- **Panel (b) added: cumulative egg-to-fry survival.** This is the panel that
  answers Reviewer 2 — the exponential models fall from ~99% at 12 °C to ~0% at
  17 °C. On daily hazard alone they look flat, which is almost certainly why
  TDM.1 was read as having no temperature response.
- **Operational range shaded, and it is the strongest new fact here.** Computed
  from the CE-QUAL-W2 Hazel Avenue series (n = 307,440 daily values, Oct–Nov, all
  alternatives and met years): **14.1–18.2 °C (5th–95th percentile), median
  16.2 °C.**

That range **straddles both crossovers (16.6 and 17.3 °C)**. The Lower American
River in October–November routinely operates in exactly the band where the choice
of TDM model reverses which model is more pessimistic. That is a much better
argument than "the models differ" and it is worth putting in the response
explicitly.

Suggested caption addition: *"Shaded band is the 5th–95th percentile of modelled
October–November daily temperature at Hazel Avenue across all alternatives and
meteorological years (14.1–18.2 °C). Dotted lines in (b) mark where each
exponential model's cumulative survival falls below Martin et al. (2017)."*

---

## Task 9 — Projected vs observed abundance

**The gap is real, it is large, and I can name the mechanism. This is the item
most likely to be pressed by a reviewer.**

NB trajectory (model-averaged, all met years):

| Projection year | 1 | 11 | 26 | 51 | 76 | 100 |
|---|---|---|---|---|---|---|
| Adults | 17,452 | 11,709 | 9,285 | 8,090 | 7,675 | 7,536 |

- Final-20-year median: **7,602**
- Observed 2014–2024 mean: **22,491** → projections settle at **34%** of it
- Observed 2011–2024 mean: 26,399 → **29%**
- Still declining at the end of the run (−5.3 spawners/yr over the final 20 years)

### Why

**The PHABSIM K update inverted the result.** K rose from 12,493 to 33,185
(+166%), but the final-20 NB median *fell* from 17,245 to 7,602 (−56%). The
mechanism: with density dependence much weaker at the higher K, the calibration
had to drop SAR from 0.004161 to 0.002685 (−35%) to still reproduce observed
escapement over 2011–2024. The lower SAR then sets a much lower unfished
equilibrium in the forecast, where temperatures are worse than the calibration
period.

Secondary contributor: forecast temperatures are warmer than the calibration mix
(mean egg survival 0.962 → 0.924, NB / exp_WF).

### ⚠️ Added after Task 12: the gap is not uniform across TDM models

The 34%-of-observed figure is an artifact of model-averaging one plausible
projection with one that collapses. Broken out (NB, median of final 20 years,
against the observed 2014–2024 mean of 22,491):

| TDM model | Panel weight | NB adult index | % of observed |
|---|---|---|---|
| Bratovich (2020) | 0.51 | 13,178 | **59%** |
| Bartholow & Heasley (2006) | 0.24 | 3,655 | 16% |
| Martin et al. (2017) | 0.25 | **6.7** | **0.03%** |
| Model-averaged | — | 7,600 | 34% |

Under the model carrying the majority of the panel weight the projection sits at
59% of recent observed escapement — a defensible gap for a fixed-climatology
projection with no hatchery supplementation. The headline 34% is dragged down by
the Martin variant projecting functional extinction under *every* alternative
including no-bypass.

**This is much easier to defend than the pooled number**, and it should be stated
this way. It also pre-empts the obvious follow-up: the Martin projection's
collapse is not an artifact of any bypass alternative, since no-bypass collapses
too.

**Ruled out:** the 1,000 cfs reference flow is *not* the problem. K at 1,000 cfs
(33,185) is already near the WUA maximum (35,171 at 1,500 cfs), so the flow
assumption is not binding — the low calibrated SAR is.

### ⚠️ Separate and serious: the calibration fit statistics cannot be regenerated

The reported **R² = 0.72, RMSE = 8,240 spawners over 2014–2024 cannot be
reproduced from the current repo.** `precompute.R:1041-1042` reads:

```r
# Since we're not calibrating, create empty calib_pred_by_variant for compatibility
calib_pred_by_variant <- list()
```

and `calib_pred_by_variant.rds` is 0 bytes. Those statistics come from an earlier
version of the pipeline. Either restore the calibration-prediction step and
recompute them, or remove them from the manuscript — as it stands they cannot be
defended if asked.

**Also note for anyone reading `results_full.rds`:** it is a single continuous
114-year projection seeded from the observed 2022–2024 escapement, but its `year`
column is labelled 2011–2124. The year axis is projection year, offset ~14 years
from calendar. `app.R` treats years < 2025 as history, which is not what they are.

---

## Task 10 — EVPI under both weight sets

*Script:* `analysis/evpi.R` → `output/evpi.csv`

States of nature are the three TDM models with the elicited panel weights as
probabilities. EVPI = E[max] − max[E].

| Objective weights | Normalisation | Best under uncertainty | EVPI | EVPI % |
|---|---|---|---|---|
| re-derived (0.75/0.20/0.05) | pooled | PB1 | **0.0265** | **4.64%** |
| re-derived (0.75/0.20/0.05) | per-state | PB2c | 0.0126 | 1.61% |
| elicited (0.40/0.50/0.10) | per-state | PB1 | 0.0120 | 2.17% |
| elicited (0.40/0.50/0.10) | pooled | NB | 0.0000 | 0.00% |

**The published 0.026 (4.6%) reproduces exactly** — under re-derived weights with
pooled normalisation. That pins down the method that was used.

**Report the range 0.000–0.027, and note that 0.026 is the largest of the four
defensible combinations.** Under the elicited weights with the same normalisation
the EVPI is *zero*, because NB is then optimal under all three TDM models and
perfect information changes nothing.

A caveat worth stating in the response, since it is a fair criticism waiting to
happen: the **normalisation choice matters as much as the weight set.** Pooled
normalisation puts all states on one scale but lets cross-TDM-model abundance
differences swamp cross-alternative differences (so the Chinook score partly
measures which model you are in). Per-state normalisation is what the app and
Figure 5 use for ranking but rescales each state separately, which is not
strictly valid for an expected-value calculation. Neither is clean — which is
itself the best argument for reporting a range rather than a point estimate.

---

## Task 11 — References to add

*Full formatted entries:* `output/new_references.md`

All three are **required, not conditional** — Panelist 5's justification, which
Table S2-8 reproduces verbatim, cites all three by name.

- Anderson, J. J., W. N. Beer, J. A. Israel, and S. Greene. 2022. Targeting river
  operations to the critical thermal window of fish incubation: Model and case
  study on Sacramento River winter-run Chinook Salmon. *River Research and
  Applications* 38(5):895–905.
- Bloomer, J., J. J. Anderson, D. Sear, S. Greene, D. Gantner, and C. Hanson.
  2023. Gastrulation and hatch as critical thermal windows for salmonid embryo
  development. *River Research and Applications* 39(1):46–53.
- Martin, B. T., P. N. Dudley, N. S. Kashef, D. M. Stafford, W. J. Reeder,
  D. Tonina, A. M. Del Rio, J. Scott Foott, and E. M. Danner. 2020. The
  biophysical basis of thermal tolerance in fish eggs. *Proceedings of the Royal
  Society B: Biological Sciences* 287(1937):20201550.

Two extras found while reading the justifications:
- Panelist 2 also cites **Anderson et al. 2018** — check it is in the reference list.
- Panelist 5 writes "Bloomer et al 2022"; the correct year is **2023**. Keep the
  panelist's text verbatim in S2-8 and cite 2023 in the reference list, or add a
  bracketed correction.

Cover letter: Anderson (2022) and Bloomer (2023) are both in *River Research and
Applications*, and both bear directly on the reviewers' concern about the thermal
window and stage-specific egg mortality.

---

## Task 12 — Front-loading mechanism ✅ VERIFIED (all four sub-checks)

*Script:* `analysis/frontloading_cohort_decomposition.R`
*Outputs:* `output/frontloading_{crossover_dates,survival_by_spawn_date,incubation_overlap,cohort_decomposition,index_by_variant}.csv`,
`figures/frontloading_cohort_decomposition.png`

**All four sub-checks now pass. The mechanism can be asserted, not merely
proposed.** Sub-checks 2 and 3 are reported below the earlier two.

### Sub-check 4 — CONFIRMED

PB3 = 9,115.5, PB5 = 9,079.8. **Difference 35.7 adults.** The manuscript's "36"
is right. Identical 21.4 Mm³, different schedule, near-identical outcome — the
honest counterweight stands as written.

### Sub-check 1 — CONFIRMED, with two corrections to Section 3.1.1

Hazel Avenue daily temperatures averaged across the four met years, differenced
against NB. Crossover = first day the alternative becomes warmer than NB and
stays warmer through Nov 30:

| Alt | Schedule ends | Crossover | Nov 16–30 mean Δ vs NB | Max Δ |
|---|---|---|---|---|
| PB1 | Nov 14 | **Nov 14** | +0.16 °C | +0.20 |
| PB3 | Nov 14 | **Nov 14** | +0.27 | +0.35 |
| PB6 | Nov 14 | **Nov 15** | **+0.41** | **+0.50** |
| PB4 | Nov 21 | Nov 21 | +0.17 | +0.43 |
| PB5 | Nov 21 | Nov 21 | +0.13 | +0.38 |
| PB2 | Nov 30 | Nov 22 | +0.03 | +0.21 |
| PB2b | — | Nov 24 | −0.13 | +0.27 |
| PB2c | — | Nov 27 | −0.24 | +0.23 |

**Correction 1:** *all eight* bypass alternatives end up warmer than NB in late
November, not just PB1, PB3 and PB6. Section 3.1.1 understates this.

**Correction 2, and it is the supportive part:** the crossover date tracks the
schedule end-date exactly as the front-loading story predicts. The three
alternatives whose bypass ends Nov 14 cross over on Nov 14–15; those running
later cross over later; PB2c, which runs latest, crosses last. **PB6 has the
largest late-November warm anomaly of any alternative**, consistent with the
"45% more water than PB4, fewer fish" pattern.

**Magnitude caveat that must be stated:** the late-November warming is only +0.03
to +0.41 °C in half-month means, while the Oct 16 – Nov 15 cooling benefit is
−0.35 to −0.91 °C — roughly two to three times larger. Direction and ordering
support the mechanism; magnitude alone does not establish that it is what drives
the abundance ranking. (Sub-check 2 below does establish it, on the hazard scale
rather than the temperature scale.)

*Both columns above are regenerated by
`analysis/frontloading_cohort_decomposition.R` Part 1 →
`output/frontloading_crossover_dates.csv`. The Nov 16–30 column reproduces the
earlier hand computation exactly; the Oct 16 – Nov 15 range supersedes an earlier
figure of −0.26 to −0.99, which had no script behind it.*

### Sub-check 2 — CONFIRMED

Both TDM families are additive in the log (`-log S = Σ h(T_d)` over the
incubation window), so the change in log-survival relative to NB splits
*exactly* over calendar days. Splitting each alternative at its own crossover
date gives the cooling benefit earned before it and the warming penalty paid
after it:

| Alt | Crossover | Benefit before (Δ −log S) | Penalty after | **Net** | Penalty as % of benefit |
|---|---|---|---|---|---|
| PB2c | Nov 27 | −0.1496 | +0.0118 | **−0.1378** | 8% |
| PB2b | Nov 24 | −0.1306 | +0.0220 | **−0.1086** | 17% |
| PB2 | Nov 22 | −0.0896 | +0.0179 | **−0.0717** | 20% |
| PB5 | Nov 21 | −0.0674 | +0.0498 | **−0.0176** | 74% |
| PB4 | Nov 21 | −0.0641 | +0.0561 | **−0.0080** | 87% |
| PB1 | Nov 14 | −0.0244 | +0.0374 | **+0.0130** | **153%** |
| PB3 | Nov 14 | −0.0504 | +0.0645 | **+0.0141** | **128%** |
| PB6 | Nov 15 | −0.0718 | +0.0928 | **+0.0210** | **129%** |

Negative = less mortality than NB. **For PB1, PB3 and PB6 the late-November
penalty exceeds the October cooling benefit outright**, which is why these three
fall below no-bypass under Martin. The ordering tracks the schedule end-date
exactly, as front-loading predicts.

The literal question asked — what fraction of redds are still incubating at the
crossover — is **100%** for every alternative, and is uninformative: incubation
runs 60–90 days, so essentially every Oct–Nov redd is in the gravel in late
November. The hazard split above is the version of the question that has an
answer.

**The decisive validation:** this net column rank-orders the Martin-only adult
index perfectly — all 9 alternatives, in exact inverse order (Spearman ρ = −1.00,
more hazard → fewer fish):

| | PB2c | PB2b | PB2 | PB5 | PB4 | NB | PB1 | PB3 | PB6 |
|---|---|---|---|---|---|---|---|---|---|
| Net Δ −log S | −0.138 | −0.109 | −0.072 | −0.018 | −0.008 | 0 | +0.013 | +0.014 | +0.021 |
| Martin adult index | 84.8 | 49.6 | 21.4 | 10.2 | 8.0 | 6.7 | 4.5 | 4.0 | 3.6 |

The incubation-window hazard balance *is* the ranking. Nothing else needs to be
invoked.

### Sub-check 3 — CONFIRMED, with one refinement

Egg survival is the channel: PB6/PB4 ratio is **0.977 for egg survival** but
**1.0008 for pre-spawn survival**, so the pre-spawn pathway contributes nothing.
The alternatives share one redd distribution, so the difference in mean
egg-to-fry survival decomposes exactly by spawn cohort.

> **⚠️ SUPERSEDED 2026-08-18. The table below is the survival channel only — do
> not quote it.** The G1 correction removed the shared-redd-distribution premise,
> so the gap splits into a survival channel *and* a composition channel, and
> composition turns out to be 41% of it. The December cohorts, shown as 0% below,
> are in fact **47%** of the gap. Current table:
> `MANUSCRIPT_REVISION_HANDOFF.md` §2.4; numbers:
> `output/frontloading_cohort_two_channel.csv`.

**PB6 − PB4**, total −0.0076 in mean egg-to-fry survival:

| Cohort | Redd share | S (PB6) | S (PB4) | Contribution | Share of gap |
|---|---|---|---|---|---|
| Oct 1–15 | 1.4% | 0.567 | 0.566 | 0.00000 | 0% |
| Oct 16–31 | 8.6% | 0.051 | 0.053 | −0.00018 | 2% |
| **Nov 1–15** | 25.0% | 0.169 | 0.184 | **−0.00369** | **49%** |
| **Nov 16–30** | 35.5% | 0.534 | 0.544 | **−0.00368** | **49%** |
| Dec 1–15 | 23.0% | 0.959 | 0.959 | 0.00000 | 0% |
| Dec 16–Jan | 6.4% | 1.000 | 1.000 | 0.00000 | 0% |

**PB6 − NB**, total −0.0190: Nov 16–30 contributes **80%**, Nov 1–15 **26%**, and
Oct 16–31 is **−6%** — that cohort is *better off* under PB6, the cooling benefit
showing up exactly where the mechanism says it should.

**The refinement:** the deficit is emphatically *not* spread evenly — 97.5% of it
comes from November-spawning cohorts holding 60.5% of the redds, and December and
January cohorts contribute exactly zero. But "late-spawning cohorts" is the wrong
label and should not be used in the manuscript. December spawners are the latest
and are entirely unaffected, because their incubation runs in water already below
the 12.14 °C Martin threshold. **The affected group is November spawners** — those
whose incubation window covers the late-November warm anomaly. Phrase it that way.

### Reproducibility caveat on this task

`egg_summary.rds` could not be reproduced from the saved artifacts: the simulated
redd set behind it was never written to disk, and `app_data/sim_redds.rds` is a
stale artifact from an older pipeline version (recomputing through it gives
correlation 0.40 against the saved values). The draw itself is seeded and would
reproduce on a full re-run of `precompute.R`, but re-running the whole pipeline
was not an option here without invalidating the numbers already reported.
Sub-checks
2 and 3 therefore avoid the sampled redd set entirely — survival is computed
deterministically for *every* spawn date from the CE-QUAL-W2 temperatures, and the
cohort weights are applied afterwards from two independent spawn-date distributions
(observed carcass surveys 2011–2024, and the stale simulated set). **Every number
above is materially identical under both weightings** — the PB6−PB4 gap is −0.0076
vs −0.0074, and the cohort shares agree to one percentage point. See the
reproducibility section below.

> **⚠️ Narrow this claim (2026-08-18).** It establishes that the *survival*
> channel does not depend on the choice of weighting, which still holds. It does
> not establish that the decomposition is complete. Both weightings apply a
> single distribution to both alternatives, so both are blind to the fact that
> alternatives now have *different* redd distributions — and agreement between
> two estimators blind in the same way says nothing about the channel they both
> miss. That channel is 41% of the PB6−PB4 gap. See
> `MANUSCRIPT_REVISION_HANDOFF.md` §2.4.

---

## Consolidated list of numbers to change

> **⚠️ The five Section 3.2 rows below were superseded by the G1 correction on
> 2026-08-18** — they are the pre-correction values. Composite scores and the
> MCDA ordering both move again. Use `OUTSTANDING_ITEMS.md` §H, which supersedes
> this whole list wherever the two disagree. The rest of the table (Figure 3
> caption, SI, EVPI, §3.1.1) is unaffected by G1 and still stands.

| Location | Currently says | Should say |
|---|---|---|
| ~~Section 3.2~~ | ~~PB1 0.537~~ | ~~**0.573**~~ → **0.558** (G1) |
| ~~Section 3.2~~ | ~~PB2c 0.502~~ | 0.502 *(unchanged — normalisation anchor)* |
| ~~Section 3.2~~ | ~~NB 0.500~~ | 0.500 *(unchanged — normalisation anchor)* |
| ~~Section 3.2~~ | ~~PB2 0.499~~ | ~~**0.534**~~ → **0.514** (G1) |
| ~~Section 3.2~~ | ~~PB4 0.491~~ | ~~**0.530**~~ → **0.515** (G1) |
| ~~Section 3.2~~ | ~~*(order)*~~ | ~~PB1 > PB2 > PB4 > …~~ → PB1 > PB4 ≈ PB2 > PB3 > PB2c > NB > PB5 > PB2b > PB6 (G1; PB4/PB2 are tied within noise) |
| Figure 3 caption | Martin onset 12.8 °C | **12.14 °C** |
| SI S2.7 | `rear_surv` = 0.5419 (calibration estimate) | **0.543 calibrated; 0.5419 is the optimizer start** |
| SI Table S2-7 | SAR 0.0025, range 0.0018–0.0034 | label as **95% bootstrap CI on the mean**; split input vs calibrated |
| SI Table S2-7 | — | add SAR calibrated **0.00269**, `rear_surv` calibrated **0.543** |
| SI S2.5 | 0–100 scoring on four criteria + revision round | **single round of direct 0–1 allocation** |
| SI S1 | Hazel `[TBD]` | cite the **WTMP deliverable** |
| Response (EVPI) | 0.026 (4.6%) | **range 0.000–0.027**; 0.026 is the upper bound |
| Section 3.1.1 | PB1, PB3, PB6 warmer than NB in late Nov | **all eight** bypass alternatives are |
| Text | observed SAR range 0.05%–0.68% | ⚠️ **could not reproduce — verify source** |
| Discussion (front-loading) | proposed mechanism | **demonstrated** — assert it (Task 12) |
| Discussion / response | "late-spawning cohorts" | **"November-spawning cohorts"** — December spawners are unaffected |
| Response (Task 1 corner) | "no panelist proposed it" | add: **all nine alternatives are functionally extinct under 100% Martin** |
| Section 3.3 / Discussion | projections are 34% of observed | add the **per-model breakdown**; under Bratovich it is 59% |

---

## Open items requiring a decision

1. **Re-run everything from one commit, or adopt the current-run numbers
   throughout.** The mixed-run problem is the biggest exposure.
2. **The R² = 0.72 / RMSE = 8,240 calibration statistics cannot be regenerated.**
   Restore the calibration-prediction step or remove the claim.
3. **The projected equilibrium is 34% of observed recent escapement and still
   declining at year 100.** This needs an explicit paragraph, because a reviewer
   who notices will press hard. The honest framing: the projections are a
   *relative* comparison across alternatives under a fixed climatology, not an
   absolute forecast of abundance.
4. **Panelist 5's weights** — prose (45/20/35) vs recorded (40/30/30).
5. **The 0.05%–0.68% SAR range** — source unverified.
6. ~~Sub-checks 2 and 3 of the front-loading mechanism~~ — **done, mechanism
   confirmed.** The Discussion and response can assert it, subject to the wording
   change from "late-spawning" to "November-spawning" cohorts.
7. **Under Martin alone the population goes functionally extinct under every
   alternative** (3.6–84.8 spawners, still declining at year 100). This is not
   currently stated anywhere. It is a strong argument in the response but a
   reviewer will also read it as a problem with the Martin implementation — decide
   whether to lead with it or footnote it.

---

## Reproducibility problems found in `precompute.R`

Three, in descending order of seriousness. None invalidates a published number,
but the first blocks independent verification and should be fixed before the code
is released with the paper.

**1. `app_data/sim_redds.rds` was stale and misleading. FIXED.** It was not
written by `precompute.R` at all — it survived from an earlier pipeline version,
and its schema differs (it has a `section` column and site labels on the
historical rows, where the current `sim_redds` sets `site = NA` for those).
Recomputing `egg_summary` through it reproduces nothing: correlation 0.40, mean
0.547 vs 0.619. Anyone who picked it up as the model's redd set — natural, given
the name — got wrong answers. *Fix applied:* `precompute.R` now saves
`sim_redds` and `sim_future` with the other artifacts.

**2. The redd draw is seeded — an earlier note in this document said otherwise
and was wrong.** `set.seed(123)` at `precompute.R:31` covers §17, and nothing
consumes the RNG stream in parallel before it (the first `future_map` is at line
883; `plan()` at line 101 does not itself draw). The draw is therefore
reproducible by running the script start to finish in one session. What was
missing was only that it was never *saved*, so the redd set behind a given set of
results could not be recovered without a full re-run — now fixed under item 1.
The residual caveat is the ordinary one: re-running a chunk interactively shifts
the RNG stream and changes the draw.

**3. Alternative-specific spawn timing is computed and then discarded.**
`precompute.R:659` splits the redds by `mgt_alt` **and** `sim_year`, which is
correct — the CLM predicts a different spawn-date distribution under each
alternative's temperatures. Line 773 then overwrites that object with a split by
`sim_year` **only**, pooling all 36 alternatives' redds into one set per year, and
that pooled set is what `eval_year()` evaluates under every alternative. The
temperature-driven shift in spawn timing is therefore averaged away before it can
affect any result. Given lines 659–662 build the per-alternative version and
nothing else consumes it, this looks like an unintended regression rather than a
modelling choice.

*Consequence for this report:* item 3 is convenient here — a redd distribution
common to all alternatives is exactly what makes the Task 12 cohort decomposition
exact. But it is a real limitation of the published results and, if a reviewer
asks whether spawn timing responds to the bypass operations, the honest answer is
that the model computes the response and then does not use it.
