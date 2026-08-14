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

## Status

The default is deliberately the **published** behaviour: the manuscript points
at this repository, so a clone must keep reproducing the published numbers
until the paper is accepted. Flipping `ARG_G1_ALT_SPAWN` to `1` by default is
the first commit of the next power-bypass cycle.

Whether to correct the numbers during the current revision instead is an open
decision for the co-authors — the sign reversal above is the fact that bears on
it most directly.
