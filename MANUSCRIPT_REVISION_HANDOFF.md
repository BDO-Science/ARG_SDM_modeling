# Folsom Bypass revision — manuscript edit package

**Manuscript:** "Balancing salmon conservation and hydropower", major revision at
*River Research and Applications*
**Prepared:** 2026-07-28
**Analysis source of truth:** `ARG_SDM_modeling` repo at commit `26d61b8`
**Companion:** `REVISION_FINDINGS.md` (the same work organised as an
analysis record, with method detail and script names)

---

## How to use this document

This is a self-contained package for editing the manuscript, the SI and the
response to reviewers. It is organised by *where the edit goes*, not by which
analysis produced it. Everything needed to make an edit is stated inline — you
should not need the repo or the companion document.

Three things it will not do for you, flagged as **DECISION** where they appear:

- It will not choose between two internally inconsistent model runs (§1).
- It will not resolve source discrepancies that need a human with access to
  external deliverables or to the panelists (§7).
- It will not decide how prominently to state a finding that is both a strong
  argument and an exposure (§7).

Confidence labels used below:

| Label | Meaning |
|---|---|
| **VERIFIED** | Recomputed from model output; reproduces exactly |
| **CORRECTED** | A published number is wrong; the replacement is verified |
| **UNRESOLVED** | Could not be verified from available sources; needs a human |

---

## 0. Status as of 2026-07-28, after editing the 07-28 document versions

Audited against `FolsomBypass_manuscript_2026-07-28_revised.docx`,
`FolsomBypass_SupportingInformation_2026-07-28_revised.docx` and
`Response_to_Reviewers_RiverResearchandApps_v3.docx`.

**The mixed-run problem (§1) is closed in the documents.** Every run-dependent
number now traces to commit `26d61b8`. No value from a superseded run (12,493,
0.004161, 17,245, 18,396, 18,953, 0.537, 0.499, 0.491) survives anywhere in
either document. The one exception is the calibration fit statistics, which come
from a pipeline version that no longer exists — see §7 item 2.

| Item | State |
|---|---|
| §2.1 Composite scores → 0.573 / 0.534 / 0.530 | applied by author |
| §2.5 Consequence table, K, calibrated SAR, `rear_surv` | already current-run |
| §5.1 Figure 3 caption 12.8 → 12.14 | applied by author |
| §3.3 SAR sentence and bootstrap-CI note | applied by author |
| §4.2 Crossovers, operational band 14.1–18.2 °C | applied by author |
| §4.3 EVPI range in the **Discussion** | applied by author |
| §4.3 EVPI in the **Results** — was still a point estimate | **edited this pass** |
| Spearman rho 0.95 → 0.96 (recomputed, 0.9624) | **edited this pass** |
| SI Table S2-7 duplicate number → renumbered S2-9 | **edited this pass** |
| SI table labelled inputs as "Calibrated Value" | **edited this pass** |
| §7.2 R² = 0.72 / RMSE = 8,240 unreproducible | **flagged in place, needs your decision** |
| §7.3 SAR range source, §7.4 Panelist 5 weights | still open |

**Edits made in this pass are highlighted bright green** in the three `.docx`
files, to distinguish them from the yellow and cyan highlighting already in the
documents. Pre-edit copies are in `manuscript_backup/`. One caveat: the Spearman
sentence already carried a yellow highlight, which is now green.

Two further problems found during the audit and fixed:

- **SI Table S2-7 was used twice** — once for the new panelist-weighting table in
  S2.5, once for "Calibrated life-cycle model parameters" in S2.8. The latter is
  now **Table S2-9**; cross-references to it still need updating.
- **That table labelled empirical inputs as calibrated values**, directly
  contradicting the S2.8 prose that says the calibrated values are 0.00269 and
  0.543. Header and cells corrected.

---

## 1. The mixed-run problem — RESOLVED in the documents, retained for the record

**This section is now history.** It is kept because it explains why several
numbers changed and which run is authoritative. See §0 for current status.

The repo was re-run twice after the MCDA numbers in Section 3.2 were written:

| Commit | What changed | K (spawners) | Calibrated SAR | NB adult index |
|---|---|---|---|---|
| `c658df3` / `e0cd8c8` | *(the run Section 3.2 was written from)* | 12,493 | 0.004399 | 18,396 |
| `16e7a70` | ATU boundaries → SacPAS v3 (400 / 558 / 958) | 12,493 | 0.004161 | 17,245 |
| `26d61b8` | PHABSIM K source updated, precompute re-run | **33,185** | **0.002685** | **7,602** |

**The consequence-table abundances (7,600–11,073) are from the current run, but
the Section 3.2 composite scores (0.537 / 0.502 / 0.500 / 0.499 / 0.491) are from
the pre-rerun run.** The two are not mutually consistent, and a reviewer
comparing them will see it. Four of the five quoted scores match `e0cd8c8` to
three decimals, which is what identifies it as the source run.

**Resolution: the current run (`26d61b8`) was adopted throughout, and the
documents now reflect it.** Verified by string audit of all three `.docx` files —
see §0. If you later re-run everything from a different commit, most numbers here
will change and this package must be regenerated.

---

## 2. Main text edits

### 2.1 Section 3.2 — composite scores are stale — **CORRECTED**

The working assumption was that the text was right and Figure 5 could not resolve
a 0.003 spread. **That is backwards: the text is stale and the figure is fine.**

| Alt | Text currently says | Replace with |
|---|---|---|
| PB1 | 0.537 | **0.573** |
| PB2c | 0.502 | 0.502 *(unchanged)* |
| NB | 0.500 | 0.500 *(unchanged)* |
| PB2 | 0.499 | **0.534** |
| PB4 | 0.491 | **0.530** |

**Replace the stated ranking with the full current one:**

> PB1 0.573 > PB2 0.534 > PB4 0.530 > PB3 0.523 > PB2c 0.502 > NB 0.500 >
> PB2b 0.482 > PB5 0.465 > PB6 0.431

Two consequences for the surrounding prose:

1. Any sentence explaining that the top alternatives are too closely spaced to
   distinguish must go. **PB1 now leads PB2 by 0.040**, which is easily visible on
   a bar chart. The old top five really were separated by ~0.003; that is no
   longer the situation.
2. If Reviewer 1 read a different order off Figure 5, they were most likely
   reading a correctly regenerated figure against stale text. Brian's "Needs to be
   fixed" is right, and **the fix is to the text, not the figure.**

Numeric labels have been added above each bar in the Figure 5 plotting code, as
the revised caption already promises.

### 2.2 Section 3.1.1 — understates how general the late-November warming is — **CORRECTED**

Currently: PB1, PB3 and PB6 can be warmer than NB in mid- to late-November.

**All eight bypass alternatives end up warmer than NB in late November**, not
three. The crossover date — the first day the alternative becomes warmer than NB
and stays warmer through Nov 30 — tracks each schedule's end date:

| Alt | Bypass schedule ends | Crossover | Nov 16–30 mean Δ vs NB | Max Δ |
|---|---|---|---|---|
| PB1 | Nov 14 | **Nov 14** | +0.16 °C | +0.20 |
| PB3 | Nov 14 | **Nov 14** | +0.27 | +0.35 |
| PB6 | Nov 14 | **Nov 15** | **+0.41** | **+0.50** |
| PB4 | Nov 21 | Nov 21 | +0.17 | +0.43 |
| PB5 | Nov 21 | Nov 21 | +0.13 | +0.38 |
| PB2 | Nov 30 | Nov 22 | +0.03 | +0.21 |
| PB2b | — | Nov 24 | −0.13 | +0.27 |
| PB2c | — | Nov 27 | −0.24 | +0.23 |

Suggested replacement sentence:

> Every bypass alternative leaves mid- to late-November temperatures at Hazel
> Avenue above the no-bypass case, with the crossover date tracking the end of
> each bypass schedule: alternatives ending 14 November cross over on 14–15
> November, while PB2c, which operates latest, does not cross over until 27
> November. The magnitude of the late-November anomaly ranges from +0.03 °C
> (PB2) to +0.41 °C (PB6) in half-month means.

**Magnitude caveat that should be stated alongside it:** the late-November
warming is +0.03 to +0.41 °C, while the 16 October – 15 November cooling benefit
is −0.35 to −0.91 °C, roughly two to three times larger. Direction and ordering
support the front-loading mechanism; magnitude alone does not — which is why the
hazard decomposition in §2.3, rather than the raw temperature anomalies, is what
actually establishes it.

### 2.3 Discussion — the front-loading mechanism can now be asserted — **VERIFIED**

Both the revised Discussion and the response to Reviewer 2 currently assert this
mechanism on the basis of inference from Table 2 schedules. **It has now been
verified against model output and can be stated as a result rather than a
conjecture** — but with one wording change that matters.

The mechanism as verified: both TDM families are additive in the log
(`−log S = Σ h(T_d)` over the incubation window), so each alternative's change in
log-survival relative to no-bypass splits *exactly* into the cooling benefit
earned before its crossover date and the warming penalty paid after it.

| Alt | Crossover | Benefit before | Penalty after | **Net** | Penalty as % of benefit |
|---|---|---|---|---|---|
| PB2c | Nov 27 | −0.1496 | +0.0118 | **−0.1378** | 8% |
| PB2b | Nov 24 | −0.1306 | +0.0220 | **−0.1086** | 17% |
| PB2 | Nov 22 | −0.0896 | +0.0179 | **−0.0717** | 20% |
| PB5 | Nov 21 | −0.0674 | +0.0498 | **−0.0176** | 74% |
| PB4 | Nov 21 | −0.0641 | +0.0561 | **−0.0080** | 87% |
| PB1 | Nov 14 | −0.0244 | +0.0374 | **+0.0130** | **153%** |
| PB3 | Nov 14 | −0.0504 | +0.0645 | **+0.0141** | **128%** |
| PB6 | Nov 15 | −0.0718 | +0.0928 | **+0.0210** | **129%** |

Units are Δ(−log S); negative means less mortality than no-bypass. For PB1, PB3
and PB6 **the late-November penalty exceeds the October cooling benefit
outright**, which is exactly why those three fall below no-bypass under Martin.

**The strongest single sentence available here:** this net column rank-orders the
Martin-only adult population index perfectly — all 9 alternatives, in exact
inverse order (Spearman ρ = −1.00, more hazard → fewer fish):

| | PB2c | PB2b | PB2 | PB5 | PB4 | NB | PB1 | PB3 | PB6 |
|---|---|---|---|---|---|---|---|---|---|
| Net Δ −log S | −0.138 | −0.109 | −0.072 | −0.018 | −0.008 | 0 | +0.013 | +0.014 | +0.021 |
| Martin adult index | 84.8 | 49.6 | 21.4 | 10.2 | 8.0 | 6.7 | 4.5 | 4.0 | 3.6 |

The balance of hazard within the incubation window *is* the ranking; no other
mechanism needs to be invoked.

**⚠️ REQUIRED WORDING CHANGE — "late-spawning cohorts" is wrong.** Both documents
currently attribute the deficit to late-spawning cohorts. The cohort
decomposition (§2.4) shows the affected group is **November spawners**. December
and January spawners are the latest of all and are affected *not at all*, because
their incubation runs in water already below the 12.14 °C Martin threshold.
Replace "late-spawning cohorts" with "November-spawning cohorts" or "cohorts whose
incubation window covers late November" everywhere it appears.

**⚠️ SUPERSEDED 2026-08-18 — the PB3 vs PB5 comparison must be rewritten, not
retained.** This entry previously read "Retain as written… the manuscript's '36'
is right, and it remains a useful honest counterweight." That was correct under
the pooled spawn-timing behaviour, which has since been corrected (G1). Under
alternative-specific spawn timing the comparison **reverses sign**: PB5 leads
PB3 by 261 ± 36 adults, the same sign at all five seeds. Identical bypass volume
(21.4 Mm³) on different schedules now produces a difference several times the
run-to-run noise, in the direction PB5's net hazard balance already predicted
(PB3 +0.014, PB5 −0.018).

The counterweight sentence therefore has to go: schedule is not "consequential
only at specific points in the decision space". Replacement wording and the full
number set are in `output/g1_revision_claims.md`. See `G1_FINDINGS.md` for the
measurement and `analysis/g1_revision_numbers.R` to regenerate.

### 2.4 Discussion — cohort decomposition, PB6 vs PB4 — **VERIFIED**

Supporting detail for §2.3; use as much as the Discussion has room for. Egg
survival is confirmed to be the only channel that matters: the PB6/PB4 ratio is
**0.977 for egg-to-fry survival** but **1.0008 for pre-spawn survival**.

**⚠️ Premise changed 2026-08-18.** This decomposition was exact *because* all
alternatives shared one redd distribution. The G1 correction makes redd
distributions alternative-specific, so the split below is now the
**survival-response component only**; an additional composition term
`sum_c (w_c,A − w_c,B) · Sbar_c` is not included. The qualitative conclusion is
unaffected — egg survival is still the only channel that matters, and the
affected group is still November spawners — but the percentages below should be
described as approximate, or the composition term quantified before they are
quoted to two figures. See `analysis/frontloading_cohort_decomposition.R` §3.

**PB6 − PB4**, total −0.0076:

| Cohort | Redd share | S (PB6) | S (PB4) | Contribution | Share of gap |
|---|---|---|---|---|---|
| Oct 1–15 | 1.4% | 0.567 | 0.566 | 0.00000 | 0% |
| Oct 16–31 | 8.6% | 0.051 | 0.053 | −0.00018 | 2% |
| **Nov 1–15** | 25.0% | 0.169 | 0.184 | **−0.00369** | **49%** |
| **Nov 16–30** | 35.5% | 0.534 | 0.544 | **−0.00368** | **49%** |
| Dec 1–15 | 23.0% | 0.959 | 0.959 | 0.00000 | 0% |
| Dec 16–Jan | 6.4% | 1.000 | 1.000 | 0.00000 | 0% |

**PB6 − NB**, total −0.0190: Nov 16–30 contributes **80%**, Nov 1–15 **26%**, and
Oct 16–31 contributes **−6%** — that cohort is *better off* under PB6, the cooling
benefit appearing exactly where the mechanism predicts.

Suggested sentence:

> The deficit is not spread evenly across the run. Cohorts spawning in November
> hold 60.5% of redds but account for 97.5% of PB6's egg-to-fry survival deficit
> relative to PB4, while October cohorts are slightly better off under PB6 and
> December and January cohorts are unaffected — their incubation occurs in water
> already below the Martin threshold.

### 2.5 Section 3.3 / Discussion — projected vs observed abundance — **VERIFIED**

The projections settle well below recent observed escapement and a reviewer who
notices will press on it. The gap is real, but it is **not uniform across TDM
models**, and stating the breakdown is far easier to defend than the pooled figure.

NB, median of the final 20 forecast years, against the observed 2014–2024 mean of
22,491:

| TDM model | Panel weight | NB adult index | % of observed |
|---|---|---|---|
| Bratovich (2020) | 0.51 | 13,178 | **59%** |
| Bartholow & Heasley (2006) | 0.24 | 3,655 | 16% |
| Martin et al. (2017) | 0.25 | **6.7** | **0.03%** |
| Model-averaged | — | 7,600 | 34% |

Under the model carrying the majority of the panel weight the projection sits at
59% of recent observed escapement — a defensible gap for a fixed-climatology
projection with no hatchery supplementation. The headline 34% is dragged down by
the Martin variant, which projects functional extinction under **every**
alternative including no-bypass — so the collapse is not an artifact of any
bypass operation.

Mechanism for the overall decline, if asked: the PHABSIM K update raised K from
12,493 to 33,185 (+166%), which weakened density dependence, which forced the
calibration to drop SAR from 0.004161 to 0.002685 (−35%) to still reproduce
observed escapement over 2011–2024. The lower SAR sets a much lower equilibrium in
the forecast, where temperatures are worse than the calibration period. **Ruled
out:** the 1,000 cfs reference flow is not the constraint — K at 1,000 cfs (33,185)
is already near the WUA maximum (35,171 at 1,500 cfs).

**Recommended framing, which should appear explicitly:**

> These projections are a *relative* comparison across alternatives under a fixed
> climatology, not an absolute forecast of abundance.

**⚠️ Related blocker — see §7, item 2:** the reported calibration fit statistics
(R² = 0.72, RMSE = 8,240 spawners) cannot be regenerated from the current
pipeline and cannot be defended if challenged.

---

## 3. Supplementary Information edits

### 3.1 SI S1 — Hazel Avenue regression statistics — **UNRESOLVED, cite the source**

The `[TBD]` placeholders cannot be filled from the repo. There is no monthly
multiple regression for Hazel Avenue *or* Watt Avenue anywhere in the
`ARG_SDM_modeling` repo, the sibling copy at
`Documents/R/reclamation/gitlab/american_river_SDM`, or any project spreadsheet
(`ARG_LAR_TempModeling_*.xlsx`, `SDM Power Bypass Temperature Modeling
Results.xlsx` — all contain scenario temperature series only, no regression
sheets). The repo consumes NWIS gauge observations (station 11446500 Hazel,
11446980 Watt) and CE-QUAL-W2 scenario output; it never fits the regression.

The Watt Avenue R² = 0.94 / RMSE = 0.76 °C / MAE = 0.58 °C must therefore come
from the WTMP deliverable.

**Action:** replace the `[TBD]` placeholders with a citation to the WTMP
deliverable, or request the Hazel statistics from whoever produced it. **Do not
compute a substitute** — it would not be the same model, and presenting it as
though it were would be worse than the current gap.

### 3.2 SI S2.7 and S2.8 — `rear_surv` — **VERIFIED, your edit was right**

**0.543 is the calibrated value; 0.5419 is the optimizer starting value.**
`precompute.R:983-991` reads:

```r
opt_combined <- optim(
  par     = c(0.0025, 0.5419),      # Starting values for SAR and rear_surv
  fn      = combined_sse, ...
  lower   = c(0.0001, 0.01),
  upper   = c(0.1, 1.0)
)
```

`calib_results.rds` holds SAR_mean = 0.002685214 and rear_surv = 0.5427946, the
latter rounding to 0.543. Both are used for all three TDM variants.

**Action:** keep S2.7 as edited. Correct S2.8 if it disagrees.

### 3.3 SI Table S2-7 — SAR row is ambiguous — **VERIFIED (partly)**

The source is `SalmonCountR/app_data/SAR LAR Releases.xlsx`. Filtering to the
**27 AMERICAN R release groups** gives mean = 0.002501 and SD = 0.002371 —
exactly the model's `SAR_mean` 0.0025 and `SAR_sd` 0.00237.

**The 0.0018–0.0034 range is a 95% bootstrap percentile CI on the mean**
(reproduces as 0.00168–0.00341; mean ± 1.96 SE gives 0.00160–0.00339). It is an
uncertainty interval on a mean, *not* a spread of observations, so it is not
comparable to an interannual range and the table must say which it is.

**Replace the SAR row with this structure**, separating empirical inputs from
calibrated estimates:

| Parameter | Value | Basis |
|---|---|---|
| SAR, input | 0.0025 (SD 0.00237) | mean of 27 American River CWT release groups |
| SAR, 95% CI on mean | 0.0018–0.0034 | bootstrap percentile CI |
| SAR, calibrated | 0.00269 | optimizer result |
| `rear_surv`, input | 0.5419 | optimizer starting value |
| `rear_surv`, calibrated | 0.543 | optimizer result |

**⚠️ Two further errors in the same SI sentence — see §7, item 3.** The sentence
currently reads:

> The 2011–2024 mean SAR was 0.25% (0.0025), with substantial interannual
> variation (range: 0.05% to 0.68%).

**(a) The range 0.05%–0.68% does not come from this file.** Twelve aggregations
were tried and none reproduces it; the closest are 0.02%–0.56% (pooled brood-year
SAR) and 0.02%–0.62% (individual release groups). The stated range is shifted
upward at both ends. The mean *does* reproduce exactly, so the file is the right
source and the range was computed from something else.

**(b) The window "2011–2024" is wrong.** The CWT file covers **brood years
2008–2019** (releases 2009–2020), and the 0.0025 mean is taken over all 27
American River release groups including broods 2008–2010 — outside the stated
window entirely. There is no 2020–2024 CWT data in the file.

**Suggested replacement, all values verified:**

> The mean SAR across 27 American River coded-wire-tag release groups (brood
> years 2008–2019) was 0.25% (0.0025, SD 0.0024). Interannual variation was
> substantial: pooled brood-year SARs ranged from 0.02% (brood year 2017) to
> 0.56% (brood year 2016).

Pooled brood-year SARs, American River releases, for reference:

| Brood year | Released | Expanded returns | SAR |
|---|---|---|---|
| 2008 | 270,000 | 158 | 0.058% |
| 2009 | 274,514 | 1,503 | 0.548% |
| 2010 | 271,171 | 2,580 | 0.951% |
| 2011 | 3,492,113 | 7,543 | 0.216% |
| 2012 | 3,277,594 | 6,591 | 0.201% |
| 2015 | 2,770,112 | 1,588 | 0.057% |
| 2016 | 2,367,561 | 13,202 | 0.558% |
| 2017 | 1,336,727 | 323 | 0.024% |
| 2018 | 2,602,318 | 3,389 | 0.130% |
| 2019 | 2,594,954 | 7,529 | 0.290% |

**Minor bug worth knowing:** in `SAR LAR Releases.xlsx` the `sar_percent` column
is identical to the `sar` column — it was never multiplied by 100. Anyone reading
that column as a percentage gets values 100× too small.

### 3.4 SI S2.5 — the elicitation protocol described is not the one used — **CORRECTED**

S2.5 currently describes panelists scoring each model 0–100 against four criteria,
normalising, and then revising after anonymous feedback. **None of that matches
the available record.** There is no criteria-scoring sheet and no Round 2 sheet
anywhere in the manuscript project folder or the repo — only
`ScoreSheet_ARG_TDM_Round1.xlsx`, which contains direct 0–1 allocations.

**Replace the protocol description with:** a single round of direct 0–1 weight
allocation across the three TDM models, with written justification. **Drop the
revision step entirely.**

### 3.5 SI Table S2-7 (new) — individual panelist weightings — **VERIFIED**

Anonymised as agreed. All five allocations sum to 1; column means confirmed.

| Panelist | Martin (2017) | Bratovich (2020) | Bartholow & Heasley (2006) | Sum |
|---|---|---|---|---|
| Panelist 1 | 0.25 | 0.25 | 0.50 | 1.00 |
| Panelist 2 | 0.10 | 0.90 | 0.00 | 1.00 |
| Panelist 3 | 0.20 | 0.70 | 0.10 | 1.00 |
| Panelist 4 | 0.30 | 0.40 | 0.30 | 1.00 |
| Panelist 5 | 0.40 | 0.30 | 0.30 | 1.00 |
| **Panel mean (weights used)** | **0.25** | **0.51** | **0.24** | 1.00 |

Machine-readable copy: `output/table_S2-7_panelist_weights.csv`.

### 3.6 SI Table S2-8 (new) — verbatim justifications — **VERIFIED**

Text is in `output/table_S2-8_justifications.md`, reproduced verbatim from the
scoresheet.

**⚠️ UNRESOLVED — see §7, item 4.** Panelist 5's prose and recorded numbers
disagree: the justification states Martin "received the highest weight (45%)",
Bratovich "20%" and Bartholow "35%", but the recorded row is 0.40 / 0.30 / 0.30.
The discrepancy becomes visible the moment the justification is reproduced
verbatim next to the table.

Also note: Panelist 5 writes "Bloomer et al 2022"; the correct year is **2023**.
Keep the panelist's text verbatim and cite 2023 in the reference list, or add a
bracketed editorial correction.

---

## 4. Response to reviewers

### 4.1 Reviewer 2, TDM weight sensitivity — the headline result — **VERIFIED**

Composite scores under five TDM weightings, with objective weights held at the
elicited 0.40 salmon / 0.50 hydropower / 0.10 steelhead and meteorological
scenarios at 0.25 each. **Top-ranked alternative in each column is bolded.**

| Alt | Bratovich only (1,0,0) | Bartholow only (0,1,0) | Martin only (0,0,1) | Equal (⅓,⅓,⅓) | Elicited (.51,.24,.25) |
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

**The direct answer to "does the top-ranked alternative change under 100%
Martin": yes.** It flips from PB1 to NB (0.5153), with PB2c second at 0.5017,
margin 0.0135. The composite ranking is *not* unconditionally invariant, and the
response must not claim that it is.

**The honest and still-strong version:**

> Holding Bratovich : Bartholow at the elicited 0.51 : 0.24 and sweeping the
> Martin weight from 0 to 1, PB1 remains top-ranked for Martin weight 0 through
> **0.988**. PB2c takes the top rank over 0.990–0.998, and no-bypass only at
> 1.000 exactly. The elicited Martin weight is 0.25, where PB1 leads the
> runner-up by 0.0395.

**Two arguments about the degenerate corner, and the second is much stronger than
the first:**

1. No panelist proposed a weighting anywhere near it — the lowest individual
   Bratovich + Bartholow allocation was 0.60 (Panelist 5).
2. **At 100% Martin weight every alternative is functionally extinct.** The
   entire nine-alternative range of the adult population index is 3.6 to 84.8
   spawners, and the trajectory never equilibrates — no-bypass declines
   monotonically across the whole 114-year run to about 7 spawners. The Chinook
   objective at that corner is a min-max normalisation across a range of a few
   dozen fish, so it carries almost no information: no-bypass scores 0.038 on the
   normalised Chinook axis and wins the composite purely on hydropower, where it
   scores 1.000 by construction.

   **So the corner where the preferred alternative flips is not a corner where
   no-bypass is better for salmon — it is a corner where the salmon axis has
   collapsed and hydropower decides.** This is worth stating plainly. See §7
   item 5 for the counter-consideration.

Supporting figure: `figures/tdm_weight_sensitivity.png` (panel a, composite score
by alternative under each weighting; panel b, the continuous Martin-weight sweep).

### 4.2 Reviewer 2, "TDM.1 has essentially no temperature response" — **VERIFIED**

The egg-to-fry survival table in the response reproduces from the model
implementation, not just from the SI parameter list. Recomputed by calling
`functions.R::tdm_exp()` and `tdm_lin_martin()` directly with ATU-paced stage
boundaries (hatch 400, emergence 958):

| T (°C) | Bratovich | Bartholow | Martin |
|---|---|---|---|
| 13 | 95.0 | 94.9 | 19.1 |
| 14 | 84.7 | 81.8 | 3.6 |
| 15 | 59.0 | 46.2 | 0.9 |
| 16 | 18.0 | 5.0 | 0.2 |
| 17 | 0.4 | 0.0 | 0.1 |

Maximum discrepancy against the table as drafted: **0.24 percentage points**. The
table can stand as written.

**Crossovers confirmed:** Bratovich falls below Martin at **17.25 °C**, Bartholow
at **16.59 °C** (the response says ~17.2 and ~16.6). This is the load-bearing
refutation and it holds.

**One refinement for precision:** there is also a *lower* crossover at 12.15 °C,
below which Martin predicts 100% survival while the exponential models are already
below it. The exact statement is that **Martin is harsher than the exponential
models between roughly 12.2 °C and 16.6 / 17.3 °C, and gentler outside that band.**

**The strongest new fact available for this response, and it is not yet in it:**
the Lower American River operates squarely inside that band. Computed from the
CE-QUAL-W2 Hazel Avenue series (n = 307,440 daily values, Oct–Nov, all
alternatives and meteorological years), October–November temperature is
**14.1–18.2 °C (5th–95th percentile), median 16.2 °C** — which straddles both
crossovers. Suggested addition:

> In October and November the Lower American River operates in exactly the
> temperature band where the choice of TDM model reverses which model is more
> pessimistic. The disagreement between TDM.1 and TDM.3 is therefore not a
> question of one model being insensitive to temperature, but of which model is
> harsher within the operationally relevant range.

### 4.3 Reviewer 1, EVPI rationale — report a range — **VERIFIED**

The published 0.026 (4.6%) reproduces exactly, under re-derived objective weights
with pooled normalisation — which pins down the method that was used. But it is
the largest of four defensible combinations:

| Objective weights | Normalisation | Best under uncertainty | EVPI | EVPI % |
|---|---|---|---|---|
| re-derived (0.75/0.20/0.05) | pooled | PB1 | **0.0265** | **4.64%** |
| re-derived (0.75/0.20/0.05) | per-state | PB2c | 0.0126 | 1.61% |
| elicited (0.40/0.50/0.10) | per-state | PB1 | 0.0120 | 2.17% |
| elicited (0.40/0.50/0.10) | pooled | NB | 0.0000 | 0.00% |

**Report the range 0.000–0.027 and note that 0.026 is the upper bound.** Under the
elicited weights with pooled normalisation the EVPI is *zero*, because no-bypass
is then optimal under all three TDM models and perfect information changes nothing.

**A caveat worth stating pre-emptively, since it is a fair criticism waiting to
happen:** the normalisation choice matters as much as the weight set. Pooled
normalisation puts all states on one scale but lets cross-TDM-model abundance
differences swamp cross-alternative differences, so the Chinook score partly
measures which model you are in. Per-state normalisation is what the app and
Figure 5 use for ranking, but it rescales each state separately, which is not
strictly valid for an expected-value calculation. Neither is clean — which is
itself the best argument for reporting a range rather than a point estimate.

---

## 5. Figures

### 5.1 Figure 3 — replotted — **VERIFIED**

- **The plotting code already used 12.14 °C**, matching
  `functions.R::tdm_lin_martin`. The caption's 12.8 was the error; the correction
  to 12.14 is right.
- x range extended to **10–18 °C**.
- **Panel (b) added: cumulative egg-to-fry survival.** This is the panel that
  answers Reviewer 2 — the exponential models fall from ~99% at 12 °C to ~0% at
  17 °C. On daily hazard alone they look flat, which is almost certainly why
  TDM.1 was read as having no temperature response.
- Operational range shaded (see §4.2).

**Suggested caption addition:**

> Shaded band is the 5th–95th percentile of modelled October–November daily
> temperature at Hazel Avenue across all alternatives and meteorological years
> (14.1–18.2 °C). Dotted lines in (b) mark where each exponential model's
> cumulative survival falls below Martin et al. (2017).

File: `figures/figure3_tdm_curves.png`; data: `output/figure3_etf_survival_by_temp.csv`.

### 5.2 Figure 5 — numeric bar labels added — **DONE**

`figures/mcda_composite_scores.png`, regenerated with labels above each bar as
the revised caption promises. Values: `output/mcda_composite_scores_v2.csv`.

### 5.3 New figure — TDM weight sensitivity — **NEW, for the response or SI**

`figures/tdm_weight_sensitivity.png`. Panel (a) composite score by alternative
under each of the five TDM weightings, top-ranked alternative highlighted; panel
(b) the continuous Martin-weight sweep showing PB1 holding the top rank to 0.988.

### 5.4 New figure — front-loading mechanism — **NEW, suggest SI**

`figures/frontloading_cohort_decomposition.png`, five panels: (a) Hazel Avenue
daily temperature minus no-bypass; (b) Martin egg-to-fry survival by spawn date
under no-bypass; (c) the same curves minus no-bypass, showing November cohorts
losing and October cohorts gaining; (d) the hazard split at each crossover date;
(e) the PB6−PB4 cohort decomposition.

Panel (c) is the one that carries the argument visually — it shows the sign
reversal by spawn date directly, and that December onward is exactly zero.

---

## 6. References to add

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

Two more found while reading the justifications:

- Panelist 2 also cites **Anderson et al. 2018** — confirm it is in the reference
  list.
- Panelist 5's "Bloomer et al 2022" should be **2023** (see §3.6).

**Cover letter:** Anderson (2022) and Bloomer (2023) are both in *River Research
and Applications*, and both bear directly on the reviewers' concern about the
thermal window and stage-specific egg mortality. Worth noting.

---

## 7. Open decisions — these need a human

Listed in the order they block work.

**1. Which model run is authoritative.** See §1. Everything else assumes the
current run (`26d61b8`).

**2. The calibration fit statistics — NOW REGENERATED, and they are much worse
than published.**

*Script:* `analysis/calibration_fit_statistics.R` → `output/calibration_fit_statistics.csv`,
`output/calibration_predictions.csv`, `figures/calibration_observed_vs_predicted.png`

The calibration-prediction step that `precompute.R` §33 stopped producing is
restored, and `app_data/calib_pred_by_variant.rds` is no longer 0 bytes. The
statistics are now reproducible from saved artifacts without re-running
`precompute.R`, so no published number was disturbed.

**The reproduction is exact:** re-running the optimiser on the reconstructed
objective returns the saved parameters to machine precision (SAR 0.0026852,
`rear_surv` 0.5427946; difference 0.00e+00 on both).

| Statistic, TDM-weighted, 2014–2024 | SI says | Regenerated |
|---|---|---|
| R² (squared correlation) | 0.72 | **0.13** |
| R² (Nash–Sutcliffe) | — | **−0.59** |
| RMSE | 8,240 | **13,345** |
| MAPE | 24% | **51%** |

**The published values cannot be recovered under any variant tested** — not over
all 14 years (R² = 0.31, RMSE = 11,829), not per-variant, and not by reverting
the carrying capacity to the pre-rerun 12,493 (R² = 0.064, RMSE = 11,812,
MAPE = 53%). They belong to a structurally different earlier pipeline.

A negative Nash–Sutcliffe means the model predicts observed escapement worse than
the mean of the observations would. The model produces a nearly flat
17,000–25,000 trajectory and misses the 2023–2024 surge entirely (2024: observed
45,541, predicted 16,942).

**This also falsifies a claim in the same SI paragraph** — that the model
"captured the magnitude and general temporal dynamics… including major peaks in
2015, 2018, and 2023–2024." 2023 and 2024 are the two worst-fit years in the
record.

**Recommendation: replace the statistics with the regenerated values and reframe
the claim** — the model is a *relative* comparison across alternatives under a
fixed climatology, not a predictive hindcast; two parameters are fitted to
reproduce general magnitude, and interannual variation is driven by ocean
conditions the model does not represent. The alternative is to delete the fit
statistics and the temporal-dynamics claim altogether and describe the
calibration as magnitude-matching. **Do not quote R² = 0.72 — it is not
reproducible and the code is public.**

**3. The 0.05%–0.68% observed SAR range.** Not reproducible from
`SAR LAR Releases.xlsx` under any of twelve aggregations, and the accompanying
window "2011–2024" does not match the data, which ends at brood year 2019 (see
§3.3). The mean 0.25% reproduces exactly, so the file is the right source. Either
find where the range came from, or adopt the verified replacement sentence in
§3.3. **Decide which** — the replacement is safe, but if the original range came
from a broader dataset you may prefer to cite that instead.

**4. Panelist 5's weights.** Prose says 45 / 20 / 35; the recorded row is
40 / 30 / 30. If the prose values were used the panel mean would be
**0.26 / 0.49 / 0.25** instead of 0.25 / 0.51 / 0.24. The effect on results is
small, but the discrepancy is visible in Table S2-8. Decide which is
authoritative and either correct the sheet or add an editorial note.

**5. How prominently to state the Martin collapse.** That every alternative goes
functionally extinct under 100% Martin weight (§4.1) is the strongest available
argument about the degenerate corner. It is also an invitation for a reviewer to
ask whether the Martin implementation is right, given it projects extinction under
conditions the population has historically persisted in. Both readings are
available; leading with it or footnoting it is a judgement call.

**6. Whether to disclose the code issues in the release.** Two were found in
`precompute.R`. Neither invalidates a published number. **Both are now fixed or
documented in the code itself** (commit pending).

- **`app_data/sim_redds.rds` was a stale artifact** — not written by
  `precompute.R` at all, left over from an older pipeline version with a
  different schema. Recomputing `egg_summary` through it reproduces nothing
  (correlation 0.40). Anyone picking it up as the model's redd set — natural,
  given the name — would get wrong answers. **Fixed:** `precompute.R` now saves
  `sim_redds` and `sim_future` alongside the other artifacts, so the file stays
  current and the redd set behind a given result is recoverable without a full
  re-run.

  *Note on reproducibility, corrected:* the redd draw itself **is** seeded —
  `set.seed(123)` at `precompute.R:31` covers it, and nothing consumes the RNG
  stream in parallel beforehand. The draw is reproducible by running the script
  start to finish in one session. It was simply never saved. This is a much
  milder issue than not being reproducible at all, and needs no disclosure
  beyond the usual note that the script should be run end-to-end.

- **Alternative-specific spawn timing is computed and then discarded.**
  `precompute.R:659` splits redds by alternative and year, correctly; line 773
  overwrites that with a split by year only, pooling all 36 alternatives.
  `eval_year()` then evaluates the pooled set under every alternative's
  temperatures, so the CLM's temperature-driven shift in spawn timing never
  affects any result. Nothing else consumes the line-659 object, which suggests
  it was unintended. **Left as-is and flagged with a comment in the code** —
  correcting it would change every downstream number, so it is your call, not a
  silent fix.

  **This one may warrant disclosure.** If a reviewer asks whether spawn timing
  responds to the bypass operations, the honest answer is that the model computes
  the response and does not propagate it. The defensible framing is that spawn
  timing is held common across alternatives so that differences in the results
  are attributable to incubation thermal exposure alone — which is true, and is
  also what makes the Task 12 cohort decomposition exact.

---

## 8. Analysis provenance

Every number in this document is reproducible from the repo at commit `26d61b8`
by running these scripts, except where marked UNRESOLVED:

| Script | Produces |
|---|---|
| `analysis/tdm_weight_sensitivity.R` | §4.1 |
| `analysis/mcda.R` | §2.1, §5.2 |
| `analysis/figure3_tdm_curves.R` | §4.2, §5.1 |
| `analysis/elicitation_tables.R` | §3.5, §3.6 |
| `analysis/evpi.R` | §4.3 |
| `analysis/frontloading_cohort_decomposition.R` | §2.2, §2.3, §2.4, §2.5, §4.1 (item 2), §5.4 |

**One methodological note on §2.3–2.5.** The redd set behind the published
`egg_summary.rds` was never saved (§7 item 6), so those results do not use it.
Survival is computed deterministically for *every* spawn date from the CE-QUAL-W2
temperatures, and cohort weights are applied afterwards from two independent
spawn-date distributions: observed carcass surveys 2011–2024, and the stale
simulated set. **Every result is materially identical under both** — the PB6−PB4
gap is −0.0076 vs −0.0074 and cohort shares agree to one percentage point — so
none of the conclusions depend on the choice of weighting. This is stronger
evidence than the pipeline's own single stochastic realisation would have been.
