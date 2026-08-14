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

### The top choice does not change

**PB1 ranks first in both arms at every seed.** The headline recommendation is
unaffected.

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

### The PB3/PB5 counterweight reverses sign

`REVISION_FINDINGS.md` sub-check 4 offers this as an honest counterweight in
the manuscript: PB3 and PB5 use identical 21.4 Mm³ on different schedules and
produce near-identical outcomes, "difference 36 adults."

| Arm | PB3 − PB5 | Spread across seeds |
|---|---|---|
| Legacy (published) | **+43 adults** (PB3 ahead) | ±4 |
| Fixed | **−256 adults** (PB5 ahead) | ±40 |

**The sign reverses at every seed, by roughly six times the noise.** The
corrected model does not say the two are near-identical; it says PB5 is better,
by seven times the quoted magnitude, in the opposite direction.

This moves the model *toward* the mechanism the paper argues for. PB3's bypass
ends Nov 14 with a net hazard balance of +0.014; PB5 ends Nov 21 at −0.018.
Under pooling those two produce near-identical abundance despite clearly
different hazard balances — a quiet incoherence. Under the fix, PB5 pulls ahead
in the direction its hazard balance predicts.

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
