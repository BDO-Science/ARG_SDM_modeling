# G1 — alternative-specific spawn timing: what the defect is and what it costs

Measured 2026-08-14 on branch `g1-spawn-timing`. Reproduce with
`analysis/compare_g1.R` (one pair) and `analysis/compare_g1_seeds.R`
(multi-seed).

---

## The defect

`precompute.R` fits an ordinal regression (CLM) predicting *when* fish spawn
from October/November water temperature, then simulates redds separately for
every alternative — `sim_future` is a 36-alternative × 100-year grid, each
alternative drawing from its own CLM probabilities. Temperature → spawn timing
→ which days eggs sit in the gravel → thermal mortality is a real causal
channel in the model.

That work was then discarded. On `main`, `precompute.R:659` built the
alternative-specific split and `precompute.R:773` overwrote the same variable
with a split by year only, pooling all 36 alternatives' redds. `eval_year()`
evaluated every alternative against that pooled set, so each alternative's egg
survival used a spawn-date distribution averaged across all 36 rather than its
own. Nothing read the line-659 object in between: it was built and dropped.

**Consequence:** alternatives could differ only through incubation exposure,
never through when the eggs were laid. That is a coherent controlled comparison
— differences attributable to incubation temperature alone — but it was
undocumented, and the dead code made it look accidental rather than chosen.

---

## Method

Both behaviours run from one source tree: `ARG_G1_ALT_SPAWN` selects the arm
(`0` = published pooling, the default; `1` = alternative-specific) and
`ARG_SEED` overrides `set.seed(123)`.

The fix reduces the redd sample behind each survival estimate from ~36 × N_redd
to N_redd, so the corrected pipeline is inherently noisier. A raw
legacy-vs-fixed difference therefore confounds signal with sampling noise. Both
arms were run at **5 seeds** (123, 456, 789, 1011, 1213), giving:

* **noise** — spread *within* an arm across seeds
* **effect** — paired (fixed − legacy) difference at the *same* seed

**Reproducibility control:** a fresh legacy run reproduced the committed
baseline exactly — zero difference across all 16,416 `egg_summary` rows and
identical composite scores. Every difference below is attributable to G1.

---

## Results

### The selected alternative is unaffected

**PB4 is the manuscript's selected alternative** — the decision-makers chose it
"despite it not being the top-ranked option in the MCDA", on their judgment
about balancing hydropower against fish benefit. PB1 tops the composite score
(0.573) but is not the selection, so composite rank is the wrong quantity to
check for reassurance here.

On the right quantity the result is stronger than a rank comparison would be:
**PB4's abundance is the least affected of all nine alternatives** — a mean
effect of +3.7 fish across seeds, essentially zero. And because the selection
rested on decision-maker judgment rather than MCDA rank, its basis is untouched
regardless.

For completeness, PB1 does rank first in both arms at every seed.

### The effect is invisible in abundance and unmistakable in the composite

| | Within-arm noise (SD) | G1 effect | Exceeds 2× noise |
|---|---|---|---|
| Chinook spawner metric | fixed ±173–217 fish, legacy ±15–21 | 4–318 fish | **0 of 9** |
| MCDA composite | fixed ±0.002–0.004 | 0.015–0.029 | **9 of 9**, 3.5–10.9× |

Seed noise is largely a *common* shift across all nine alternatives, and the
composite's min–max normalisation cancels common shifts. G1 changes
alternatives' standing *relative to each other*, which survives normalisation.
So the effect is real and reproducible in relative terms — and relative
standing is what the MCDA consumes — while being undetectable in absolute
abundance.

> **Caveat:** NB and PB2c show exactly zero composite change with zero noise
> *by construction*, not because they are insensitive. They are the min–max
> anchors (NB the Chinook and hydropower minimum, PB2c the Chinook maximum), so
> they are pinned at 0 and 1. The composite understates sensitivity at the
> extremes.

### One ranking change is real; two apparent ones are noise

Comparing whole rankings is too coarse — a single noisy pair makes the entire
ranking look unstable. Ordering each of the 36 pairs separately, in each arm,
at every seed:

| Pair | Legacy | Fixed | Reading |
|---|---|---|---|
| **PB2b vs PB5** | PB2b > PB5 at all 5 seeds | PB5 > PB2b at all 5 seeds | **Real G1 effect** |
| PB2 vs PB4 | PB2 > PB4 at all 5 seeds | seed-dependent | Not determinable from one run |
| PB2c vs PB3 | PB3 > PB2c at all 5 seeds | seed-dependent | Not determinable from one run |

The remaining 33 pairs are stable and unchanged.

### The PB3/PB5 counterweight reverses sign — in the paper's favour

The manuscript's Discussion currently reads:

> "...the two alternatives bypassing identical volumes on different schedules
> (PB3 and PB5) differed by fewer than 40 adults, so release schedule is
> consequential at specific points in the decision space rather than throughout
> it."

That sentence is a hedge *against the paper's own front-loading argument* — it
concedes that schedule does not always matter.

| Arm | PB3 − PB5 | Spread across seeds |
|---|---|---|
| Legacy (published) | **+43 adults** (PB3 ahead) | ±4 |
| Fixed | **−256 adults** (PB5 ahead) | ±40 |

**The sign reverses at every seed, by roughly six times the noise.**

So the correction *removes* the caveat rather than creating one. It also
resolves a quiet incoherence in the current results: PB3's bypass ends Nov 14
with a net hazard balance of +0.014, PB5 ends Nov 21 at −0.018, yet under
pooling the two produce near-identical abundance despite clearly different
hazard balances. Under the fix PB5 pulls ahead, in the direction its hazard
balance predicts. The schedule-matters mechanism gets *more* internally
consistent, not less.

### Manuscript claims that would need regenerating

| Claim | Where | Why it moves |
|---|---|---|
| PB3 vs PB5 "fewer than 40 adults"; schedule consequential only "at specific points" | Discussion | Reverses sign (above) |
| PB6 vs PB4 adult index, "9,350 against 9,545" | Discussion | Both values shift |
| Efficiency "roughly 47 additional adults per million m³ ... against 76 for PB4 and PB2c" | Discussion | Derived from the above |
| Volume vs salmon benefit, "Spearman rho = 0.95" | Discussion | Recomputed across all nine |

`MANUSCRIPT_REVISION_HANDOFF.md` §2.3 currently instructs "**Retain as
written:** the PB3 vs PB5 comparison." That instruction is superseded if the
correction is applied.

### Pooling was hiding real uncertainty

Per-alternative estimates carry ±173–217 fish of Monte Carlo noise across
seeds. Pooling suppressed that to ±15–21 and made single-run point estimates
look far more precise than they are.

**Any claim resting on a difference smaller than roughly 400 fish is not
supportable from a single run** — including in the published numbers. Report a
mean over several seeds with the spread stated, rather than one seed to one
decimal place.

---

## Status — RESOLVED, correction applied

**Decision (2026-08-18): correct during the current revision.** The
alternative-specific behaviour is now the **default**; `ARG_G1_ALT_SPAWN=0`
reproduces the superseded pooled numbers for comparison. The sign reversal above
was the deciding fact — the correction removes a hedge against the paper's own
argument rather than creating a problem, so there was no case for deferring it
to the next power-bypass cycle.

Replacement values for every affected claim are tabulated below and written to
`output/g1_revision_{numbers,contrasts,efficiency}.csv` by
`analysis/g1_revision_numbers.R`, from the five-seed snapshots. They are means
across seeds with the run-to-run range, not single-run point estimates.

> **These are measurement values, not the reported ones.** The reporting
> convention (`OUTSTANDING_ITEMS.md` §H0) is that **levels** in the paper come
> from the committed run — generate them with `analysis/reporting_values.R` — with
> the ±173–217 spread quoted alongside, while **differences** between
> alternatives are reported as the paired within-seed contrasts below. So the
> PB3−PB5 figure of −261 ± 36 *is* what the paper reports; the per-alternative
> levels in the next table are not. Subtracting two levels out of Table 3
> discards the cancellation that makes the contrasts tight and would understate
> what the model resolves.

| Claim | Published | Corrected |
|---|---|---|
| PB3 − PB5 | +43 ± 4 (PB3 ahead) | **−261 ± 36 (PB5 ahead)** |
| PB6 vs PB4 adult index | 9,352 vs 9,546 | **9,199 vs 9,550** |
| PB6 − PB4 gap | −194 ± 3 | **−352 ± 21** |
| PB2b − PB5 composite | +0.0178 (PB2b ahead) | **−0.0120 (PB5 ahead)** |
| Efficiency, PB6 / PB4 / PB2c | 47 / 76 / 76 | **45 / 79 / 84** |
| Spearman ρ, volume vs benefit | 0.9624 | **0.9205** |

PB4 — the decision-makers' selection — moves by **+3.7 adults**, the smallest
change of all nine. PB1 still ranks first in the MCDA at every seed in both
arms.

Bypass volumes were transcribed from manuscript Table 2 and are now hard-coded
in `analysis/g1_revision_numbers.R`; they were not recorded anywhere in the
repository before. Override with `G1_VOLUMES` if Table 2 changes.

**Validation.** Run with `ARG_G1_ALT_SPAWN=0`, the pipeline reproduces published
Table 3 across all nine alternatives (7,600 / 8,560 / 10,505 / 10,974 / 11,073 /
9,116 / 9,545 / 9,080 / 9,350) and the published Spearman ρ of 0.9624 to four
decimal places. The corrected values above are therefore differences from the
paper's own numbers, not from a re-derivation of them.

**Two Discussion sentences need restructuring, not just renumbering:**

1. The efficiency sentence pairs "76 for PB4 and PB2c" — those two no longer
   share a value (79 and 84).
2. A stronger claim is now available in its place: PB3 and PB5 bypass the
   identical 21.4 Mm³ and their efficiencies diverge, 72 against 84 adults per
   million m³. That makes the schedule-matters point with volume held exactly
   constant, which is what the PB3/PB5 pairing was chosen for.

The Spearman ρ has no seed spread — it is rank-based, and run-to-run noise
shifts all nine alternatives together without reordering them. That claim is
stable from a single run; the abundance claims are not.

---

## A second mechanism the correction exposes

The fix broke the premise behind the cohort decomposition, and repairing it
turned up something worth reporting in its own right.

The PB6/PB4 egg-survival gap was split by spawn cohort on the assumption that
both alternatives share one redd distribution. They no longer do, so the gap
splits exactly into two channels (`frontloading_cohort_decomposition.R` Part 3b,
identity asserted in code):

| | PB6 − PB4 | PB6 − NB | PB4 − NB |
|---|---|---|---|
| Survival response (as previously reported) | −0.01401 | −0.01341 | +0.00052 |
| **Composition (spawn-date shift)** | **−0.00983** | **−0.00510** | **+0.00482** |
| Total | −0.02385 | −0.01850 | +0.00534 |
| Composition share | **41%** | **28%** | **90%** |

The survival channel reproduces the existing Part 3 result under the simulated
weighting (−0.01392 and +0.00035), which validates the implementation. What was
missing is the composition channel, and it is not small.

**The December cohorts are the headline.** Their survival contribution is exactly
zero in every comparison — December incubation runs below the 12.14 °C Martin
threshold, so survival is identical under every alternative. But their
*contribution to the gap is large and entirely compositional*: Dec 1–15 alone
accounts for −0.00760, about a third of the PB6−PB4 gap. PB6 holds 29.3% of its
redds in the December window against PB4's 30.5%; front-loading shifts fish out
of the safe window and into risky late November.

So front-loading now damages the population two independent ways: eggs already in
the gravel see a worse thermal window, *and* more eggs get laid into that window
in the first place. That is a stronger version of the paper's own argument, and
the second half of it is currently unstated.

**This invalidates a wording instruction, not just a number.**
`MANUSCRIPT_REVISION_HANDOFF.md` §2.3 flags a "REQUIRED WORDING CHANGE" —
replace "late-spawning cohorts" with "November-spawning cohorts", on the grounds
that "December and January spawners are the latest of all and are affected *not
at all*". The first half still stands; the justification no longer does.
December spawners are affected substantially, through composition rather than
survival. Describe the channel, not just the cohort.
