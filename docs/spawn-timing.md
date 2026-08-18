# Alternative-specific spawn timing

**Changed 2026-08-18.** Egg survival is now evaluated against the redds simulated
under each alternative's own temperatures. Previously redds were pooled across
all alternatives first. This page records what the model does, why it changed,
and what the change moves.

## What the model does

Spawn timing is predicted by an ordinal regression (CLM) fitted to observed redd
construction dates against mean October and November water temperature. Warmer
autumn temperatures shift the predicted distribution of spawn dates.

Each management alternative produces a distinct temperature regime, so the CLM
produces a distinct spawn-date distribution under each. `precompute.R` simulates
that grid — 36 alternative × meteorological-year combinations — as `sim_future`,
then draws redds accordingly.

Temperature therefore reaches egg-to-fry survival through two pathways:

1. **Incubation exposure** — the temperature eggs experience once in the gravel.
2. **Spawn timing** — *when* eggs are deposited relative to the autumn cooling
   curve, which determines how much of the warm period their incubation window
   overlaps.

## What was wrong

The simulated redds were pooled across all 36 alternative × year combinations
before survival was computed, and every alternative was then evaluated against
that single pooled spawn-date distribution.

The effect was to sever pathway 2. Alternatives could differ only through
incubation exposure; the spawn-timing shift the CLM exists to represent never
reached any result. The alternative-specific split was built and then overwritten
by the pooled one, with nothing reading it in between — so this was an error, not
a modelling choice, though it is a coherent controlled comparison if read as one.

## What it changes

Both pathways are now active. Effects were measured across five random seeds in
both configurations, pairing within seed.

**Relative standing changes more than absolute abundance.** Seed-to-seed
variation is largely a common shift across all nine alternatives, so it cancels
in a paired contrast and in the min–max normalisation the composite uses.

| | Pooled | Alternative-specific |
|---|---|---|
| PB3 − PB5 adult index | +43 (PB3 ahead) | **−261 (PB5 ahead)** |
| PB2b − PB5 composite | +0.018 (PB2b ahead) | **−0.012 (PB5 ahead)** |
| Composite, PB1 / PB2 / PB4 | 0.573 / 0.534 / 0.530 | 0.555 / 0.513 / 0.514 |
| Rank correlation, bypass volume vs adult index | 0.962 | 0.921 |

PB1 remains the highest-scoring alternative under both. PB4's abundance is the
least affected of the nine, changing by about 4 fish.

**The PB3/PB5 result is the most informative.** Those two bypass identical
volumes (21.4 Mm³) on different schedules, so they isolate schedule from volume.
Under pooling they came out nearly identical despite clearly different net
thermal hazard balances during incubation (+0.014 against −0.018) — an
inconsistency with no mechanistic explanation. With spawn timing active, PB5
pulls ahead in the direction its hazard balance predicts.

## Run-to-run variability

Pooling drew each survival estimate from a redd sample roughly 36× larger than
any single alternative's, which suppressed Monte Carlo noise in the reported
values. Measured across five seeds:

| | Run-to-run SD, adult index |
|---|---|
| Pooled | ±15–21 fish per alternative |
| Alternative-specific | ±173–217 fish per alternative |

The higher figure is the honest one — the lower was an artefact of pooling, not
extra precision.

**Report differences between alternatives as paired within-seed contrasts, not
as subtractions of two levels.** The common seed shift cancels in a pairing, so
contrasts are far tighter than levels: ±21–36 against ±173–217. Subtracting two
numbers out of a summary table discards that cancellation.

## Two channels in the front-loading mechanism

Because alternatives now have different redd distributions, a difference in mean
egg-to-fry survival decomposes exactly into two terms:

    gap = Σ_c wbar_c (S_c,A − S_c,B)      survival response
        + Σ_c (w_c,A − w_c,B) · Sbar_c    composition

where `w_c` is the redd share of spawn cohort `c` and `S_c` its survival. Under
pooling the composition term was zero by construction.

It is not small. For PB6 − PB4 it is **41%** of the gap. The December cohorts are
the clearest case: their survival is *identical* under every alternative, because
December incubation runs below the 12.14 °C Martin et al. (2017) threshold — yet
Dec 1–15 alone accounts for about a third of the gap, entirely through
composition. Front-loaded release advances spawn timing and moves redds out of
the thermally safe December window into late November.

So front-loading acts twice: it worsens the thermal window for eggs in the
gravel, *and* places more eggs into that window.

`analysis/frontloading_cohort_decomposition.R` Part 3b computes this, asserting
the identity rather than assuming it.

## Reproducing either configuration

| Variable | Default | Effect |
|---|---|---|
| `ARG_SPAWN_TIMING` | `alternative` | `alternative` — each alternative uses its own simulated redds. `pooled` — the superseded behaviour, retained so earlier results can be regenerated. |
| `ARG_SEED` | `123` | Overrides `set.seed(123)` for replication across seeds. |

```sh
Rscript SalmonCountR/precompute.R                      # current behaviour
ARG_SPAWN_TIMING=pooled Rscript SalmonCountR/precompute.R   # superseded
ARG_SEED=456 Rscript SalmonCountR/precompute.R              # different draw
```

The repository state prior to this change is tagged `as-submitted-2026-07-28`.

## Scripts

| Script | Purpose |
|---|---|
| `analysis/reporting_values.R` | The reported values, from the committed run — adult index, composite and rank, volume-normalised benefit, rank correlation |
| `analysis/spawn_timing_effect.R` | Multi-seed measurement of what the change moves, as paired contrasts with spread |
| `analysis/compare_spawn_timing.R` | One pooled-vs-alternative pair, for a quick check |
| `analysis/compare_spawn_timing_seeds.R` | Full multi-seed replication, separating effect from Monte Carlo noise |

The multi-seed scripts need snapshots from several `precompute.R` runs; point
`SPAWN_TIMING_SNAPROOT` at a directory of `<mode>_seed<n>/` folders. Budget about
15 minutes per run.
