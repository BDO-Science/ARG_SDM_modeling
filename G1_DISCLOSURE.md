# G1 disclosure — draft text for the revision

Three pieces of text. Drop them in and edit for voice; the numbers are from
`output/g1_revision_claims.md` and regenerate with
`analysis/g1_revision_numbers.R`.

---

## 1. Methods — replaces/extends the spawn-timing paragraph

> Spawn timing was modelled with an ordinal regression predicting the
> distribution of redd construction dates from mean October and November water
> temperature. Because each management alternative produces a distinct
> temperature regime, the redd date distribution is simulated separately under
> each alternative, and each alternative's egg survival is evaluated against its
> own simulated redds. Temperature therefore affects egg-to-fry survival through
> two pathways: the temperature experienced during incubation, and the shift in
> when eggs are deposited relative to the autumn cooling curve.

**Note for the response letter:** the phrase "and each alternative's egg survival
is evaluated against its own simulated redds" is the substance of the change. In
the originally submitted version all alternatives' simulated redds were pooled
before survival was computed, so only the first pathway was active.

---

## 2. Response to reviewers — the disclosure itself

> During revision we identified and corrected an error in our own analysis
> pipeline. The spawn-timing submodel simulated redd date distributions
> separately for each management alternative, as described, but those
> distributions were then pooled across alternatives before egg survival was
> computed. Each alternative was consequently evaluated against a spawn-date
> distribution averaged over all alternatives rather than its own. The effect was
> to suppress one of the two pathways by which bypass operations influence
> egg-to-fry survival — the shift in spawn timing itself — leaving only the
> effect of incubation temperature.
>
> We have corrected this and regenerated all results. The correction does not
> change the ranking of alternatives in the multi-criteria analysis: PB1 remains
> the highest-scoring alternative, and PB4, the alternative selected by the
> decision-making team, is the least affected of the nine (a change of
> approximately 4 adults). One pairwise ordering changes: PB5 now scores above
> PB2b in the composite, reversing their previous order.
>
> The correction does materially change one statement in our Discussion, and in
> a direction that strengthens rather than weakens the paper's argument. We
> previously reported that PB3 and PB5 — which bypass identical volumes on
> different schedules — differed by fewer than 40 adults, and concluded that
> release schedule is consequential only at specific points in the decision
> space. Under the corrected analysis PB5 exceeds PB3 by approximately 260
> adults, consistent at every random seed tested, and in the direction predicted
> by the two alternatives' net thermal hazard balance during incubation. We have
> revised that passage accordingly. The revised text is more consistent with the
> front-loading mechanism developed in the preceding section.
>
> We also take the opportunity to report uncertainty more honestly. Because the
> corrected pipeline evaluates each alternative against its own redd sample
> rather than a pooled sample roughly 36 times larger, per-alternative estimates
> carry substantially more Monte Carlo variability than the original analysis
> implied. Across five random seeds the run-to-run standard deviation of the
> adult population index is approximately 170–220 fish per alternative, against
> 15–20 under pooling. We now report the mean across five seeds with the
> run-to-run range, and we no longer draw conclusions from differences smaller
> than roughly 400 fish. This constraint applied equally to the originally
> submitted values; pooling concealed it rather than avoiding it.

---

## 3. Supporting Information — reproducibility note

> All results in this paper were generated with `SalmonCountR/precompute.R` at
> the commit tagged in the repository. Alternative-specific spawn timing is the
> default; setting the environment variable `ARG_G1_ALT_SPAWN=0` reproduces the
> pooled-spawn-timing behaviour of the originally submitted analysis, and
> `ARG_SEED` overrides the default random seed for replication across seeds.
> Reported values are means across seeds 123, 456, 789, 1011 and 1213.

---

## 4. Optional addition to §2, if you want the efficiency change stated too

> Two further Discussion values change. The volume-normalised benefit of the
> most front-loaded alternative, PB6, falls from approximately 47 to 45
> additional adults per million cubic metres bypassed, while PB4 rises to 79 and
> PB2c to 84; those two alternatives no longer share a common value. The rank
> correlation between bypass volume and Chinook benefit falls from 0.96 to 0.92,
> the corrected analysis attributing somewhat less of the biological benefit to
> volume alone and more to release timing. The clearest illustration is the
> PB3–PB5 pair: bypassing the identical 21.4 million cubic metres, they now
> return 72 and 84 additional adults per million cubic metres respectively.

---

## What is NOT covered by this draft

- **The PB6/PB4 cohort decomposition percentages** (49% / 49% for the two
  November cohorts) are now the survival-response component only, because the
  decomposition assumed a shared redd distribution. Either describe them as
  approximate or quantify the composition term first.
- **Figures and tables** need regenerating from the corrected `app_data`, and
  Table 3's abundance column changes for all nine alternatives.
- Whether to characterise this as an **erratum-style disclosure or a routine
  revision** is an editorial judgment for the co-authors and, ultimately, the
  editor. The text above treats it as a correction found during revision, which
  is what it is.

## One argument worth making explicitly to the editor

The corrected analysis is *more* internally consistent than the submitted one,
not less. Under pooling, PB3 and PB5 produced near-identical abundance despite
clearly different net thermal hazard balances (+0.014 against −0.018) — a quiet
incoherence with no mechanistic explanation. The correction resolves it: PB5
pulls ahead, in the direction its hazard balance predicts. The same is true of
the Spearman change, which shifts explanatory weight from volume to timing,
which is the paper's own thesis. It is worth saying that the correction supports
the argument rather than merely perturbing the numbers.
