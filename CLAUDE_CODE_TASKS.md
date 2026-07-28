---
editor_options: 
  markdown: 
    wrap: 72
---

# Claude Code task list, Folsom Bypass revision

Context: manuscript "Balancing salmon conservation and hydropower" is in
major revision at River Research and Applications. Two reviewers both
challenged the temperature-dependent mortality (TDM) model weighting.
The prose revisions are done and highlighted in
`FolsomBypass_manuscript_2026-07-27_revised.docx`. Everything below
needs the repo, the model output, or the plotting code. Items are
ordered by how much the response to reviewers depends on them.

------------------------------------------------------------------------

## 1. HIGHEST PRIORITY: TDM weight sensitivity figure (new Figure)

Reviewer 2 asked for a robustness-under-structural-uncertainty analysis
and a figure showing how much the panel weighting drives the preferred
alternative. We committed to both in the response.

Produce, for all nine alternatives (NB, PB1, PB2, PB2b, PB2c, PB3, PB4,
PB5, PB6):

-   Composite MCDA score under TDM weightings `(1,0,0)`, `(0,1,0)`,
    `(0,0,1)`, equal `(1/3,1/3,1/3)`, and the elicited
    `(0.51 Bratovich, 0.24 Bartholow, 0.25 Martin)`.
-   Hold objective weights at the elicited values (0.40 salmon, 0.50
    hydropower, 0.10 steelhead) and meteorological scenarios at 0.25
    each.
-   Report the range of Martin weight over which the top-ranked
    alternative does not change.

**The single number that matters most:** does the top-ranked alternative
change under `(0,0,1)`, that is 100% Martin? Section 3.1.2 already
establishes that PB2c is best on the *salmon objective* under every TDM
model, but the *composite* ranking (PB1 was top at 0.537) has not been
tested. If the composite ranking is also invariant, say so and it goes
in the opening line of the response.

## 2. Verify the composite scores quoted in Section 3.2

Text currently reports PB1 0.537, PB2c 0.502, NB 0.500, PB2 0.499, PB4
0.491. Brian flagged this as "Needs to be fixed" because Reviewer 1 read
a different order off Figure 5. My reading is that the text is right and
the figure simply cannot resolve a 0.003 spread. Confirm the five values
against model output, then add numeric labels above the bars in the
Figure 5 plotting code. The revised Figure 5 caption already promises
those labels.

## 3. Hazel Avenue regression fit statistics

SI Section S1 reports Watt Avenue as R2 = 0.94, RMSE = 0.76 C, MAE =
0.58 C, but reports nothing for Hazel Avenue, which Reviewer 1 noticed
and which matters more since Hazel is the primary compliance point for
the upper reach and is the temperature station assigned to sections NB,
W, 1a, 1b and 2 of the carcass survey. The revised SI has `[TBD]`
placeholders. Compute R2, RMSE and MAE for the Hazel Avenue monthly
multiple regression over the 2000-2021 calibration period. If these live
only in the WTMP deliverable rather than in our repo, say so and we will
cite the source instead.

## 4. rear_surv: 0.543 or 0.5419

The SI reports both. Section S2.8 says calibration yielded
`rear_surv = 0.543`; Section S2.7 says `rear_surv = 0.5419` and calls it
a calibration estimate, but 0.5419 is also listed as the optimizer
*starting* value. I edited S2.7 to say 0.543 and label 0.5419 as the
start value. Confirm from the calibration script which value the
projections actually used, and correct whichever direction is wrong.

## 5. Table S2-7 SAR row

The table gives SAR 0.0025 with a range of 0.0018 to 0.0034. The text
gives an observed interannual range of 0.05% to 0.68%, that is 0.0005 to
0.0068. These are different quantities and the table does not label
which it is reporting (bootstrap CI on the mean? interquartile range?).
Determine what the 0.0018 to 0.0034 range is, then restructure the table
to separate empirical inputs from calibrated estimates: SAR input
0.0025, SAR calibrated 0.00269, rear_surv input 0.5419, rear_surv
calibrated 0.543.

## 6. Build the elicitation appendix tables from the scoresheet

`ScoreSheet_ARG_TDM_Round1.xlsx` has five panelists, per-model weights,
and written justifications. The response to reviewers promises two new
SI tables:

-   Table S2-7: individual weightings by panelist (anonymize as Panelist
    1 to 5 unless we get consent to name them).
-   Table S2-8: written justifications verbatim.

Column means check out: Martin 0.25, Bratovich 0.51, Bartholow 0.24. Two
things to resolve. First, SI Section S2.5 describes panelists scoring
each model 0 to 100 against four criteria and then normalizing, but the
Round 1 sheet contains direct 0 to 1 allocations with no criteria
scores. Either find the criteria-scoring sheet or correct the SI
description to match what was actually done. Second, the SI describes a
revision round after anonymous feedback; if a Round 2 sheet exists, the
final weights should come from it, not Round 1. If Round 1 is all there
is, the SI protocol description needs to drop the revision step.

## 7. Verify the egg-to-fry survival comparison used in the response

The response to Reviewer 2 contains this table, which I computed from
the SI parameters rather than from the model code. Reproduce it from the
implementation and correct any discrepancy before we submit.

```         
h(T)   = alpha * exp(beta * T)            # exponential models, T in deg C
S      = exp(-sum(h(T)))                  # cumulative over stage duration
Bratovich egg    alpha 3.408e-11  beta 1.211
Bratovich alevin alpha 1.018e-10  beta 1.241
Bartholow egg    alpha 1.475e-11  beta 1.392
Bartholow alevin alpha 2.521e-12  beta 1.461
Martin           S = exp(-0.026 * sum(max(T - 12.14, 0)))
Stage durations at constant T: egg 400/T days, alevin 558/T days (ATU 400 and 958 total)
```

Egg-to-fry survival, percent, at constant temperature:

| T (C) | Bratovich | Bartholow | Martin |
|-------|-----------|-----------|--------|
| 13    | 95.0      | 94.9      | 19.2   |
| 14    | 84.8      | 81.9      | 3.7    |
| 15    | 58.8      | 46.4      | 0.9    |
| 16    | 18.0      | 5.1       | 0.2    |
| 17    | 0.4       | 0.0       | 0.1    |

Crossovers to check: Bratovich becomes harsher than Martin above roughly
17.2 C, Bartholow above roughly 16.6 C. This crossover is load-bearing
in the response, because it refutes the reviewer's claim that TDM.1 has
"essentially no temperature response."

## 8. Replot Figure 3

Current caption said the Martin onset is 12.8 C; SI says 12.14 C. I
corrected the caption to 12.14, so check which value the plotting code
uses. While there, extend the plotted x range to at least 10 to 18 C,
shade the Lower American River October to November operational range,
and consider adding a second panel showing cumulative egg-to-fry
survival rather than daily hazard. If the current figure tops out near
15 C the exponential models look flat, which is the most plausible
reason Reviewer 2 read TDM.1 as having no temperature response.

## 9. Sanity-check projected versus observed abundance

The consequence table gives model-averaged adult population indices of
7,600 (NB) to 11,073 (PB2c), taken as the median of the final 20 years
of a 100-year projection. Calibration fit statistics are R2 = 0.72 with
RMSE = 8,240 spawners over 2014-2024, which implies observed escapement
well above the projected equilibrium. Check whether the projections are
drifting to a lower equilibrium than the calibration period, and if so
why (flow assumption in the carrying-capacity lookup? the 1,000 cfs
reference flow? seed spawners from 2022-2024?). This matters because the
response leans on the population having persisted historically, and a
reviewer who notices the gap will press on it.

## 10. EVPI under both weight sets

EVPI of 0.026 (4.6%) is computed with objective weights re-derived to
reproduce the observed ranking (0.75 salmon, 0.20 hydropower, 0.05
steelhead) rather than the elicited weights (0.40, 0.50, 0.10). Reviewer
1 asked for the rationale, which I have written up, but the
re-derivation is close to circular. Compute EVPI under both weight sets
so we can report a range instead of a single figure.

## 12. Verify the front-loading mechanism (added after the schedule-versus-volume discussion)

Both the revised Discussion and the response to Reviewer 2 now assert a
mechanism that has only been inferred, not checked against model output.
The claim is: alternatives that draw the cold-water pool down early
leave mid- to late-November temperatures elevated, and under the Martin
et al. (2017) formulation this is why several bypass alternatives fall
below no-bypass.

Evidence it currently rests on: Table 2 schedules, the
volume-versus-abundance pattern below, and the existing Section 3.1.1
statement that PB1, PB3 and PB6 can be warmer than NB in mid- to
late-November.

```         
alt    volume(Mm3)  abundance   adults per Mm3 vs NB
NB      0.0          7,600        -
PB1    12.2          8,560       78.7
PB3    21.4          9,116       70.8
PB5    21.4          9,080       69.2
PB4    25.7          9,545       75.7   <- selected
PB6    37.2          9,350       47.0   <- 45% more water than PB4, fewer fish
PB2    42.2         10,505       68.8
PB2c   45.9         11,073       75.7
PB2b   49.5         10,974       68.2   <- most water, not most fish
Spearman rho(volume, abundance) = 0.95
```

To check:

1.  Pull the CE-QUAL-W2 daily temperature series at Hazel Avenue for
    each alternative under each meteorological scenario. Confirm that
    PB6 (and PB1, PB3) exceed NB during the mid- to late-November
    window, and identify the crossover date.
2.  Cross that against the modeled spawn-date distribution to get the
    fraction of redds whose incubation window overlaps the period when
    the front-loaded alternatives are warmer than NB.
3.  Under the Martin formulation specifically, decompose the adult
    population index by spawn cohort for PB6 versus PB4 to confirm that
    the deficit comes from late-spawning cohorts rather than being
    spread evenly.
4.  Confirm the PB3 versus PB5 comparison. Identical volume (21.4 Mm3),
    different schedule, and the outcomes differ by only 36 adults. If
    that holds it is a useful honest counterweight and is already stated
    in both documents.

If the decomposition does not support the mechanism, the assertion needs
downgrading to a conjecture in both the Discussion and the response
before resubmission.

## 11. Reference additions

Add to Zotero and the reference list:

-   Anderson, J.J., W.N. Beer, J.A. Israel, and S. Greene. 2022.
    Targeting river operations to the critical thermal window of fish
    incubation: Model and case study on Sacramento River winter-run
    Chinook Salmon. River Research and Applications 38(5):895-905.
-   Bloomer, J., J.J. Anderson, D. Sear, S. Greene, D. Gantner, and C.
    Hanson. 2023. Gastrulation and hatch as critical thermal windows for
    salmonid embryo development. River Research and Applications
    39(1):46-53.
-   Martin et al. 2020, Proc. R. Soc. B 287:20201550, if we cite the
    panelist justifications that reference it.

Both Anderson and Bloomer are in this journal, which is worth noting in
the cover letter.
