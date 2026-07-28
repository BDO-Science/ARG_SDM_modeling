# SalmonCountR — options for annual data updates

**For discussion with Brian**
**Prepared:** 2026-07-28
**Question on the table:** how do new temperature scenarios and new carcass/escapement
data get into the decision-support app each year, without routing every request
through one person?

---

## The short version

The request "let non-R users load new data and get results" is really **two
different jobs**, and they should not be solved the same way.

| | Job A — explore a new temperature scenario | Job B — update the model's data foundation |
|---|---|---|
| What changes | CE-QUAL-W2 scenario temperatures | Carcass survey, escapement, CWT returns |
| What it affects | The forecast only | Calibration, spawn timing, the whole model |
| Cost to compute | **Under a second** (measured) | Minutes, plus re-fitting two models |
| Judgment required | None — it is a what-if | Real — someone must decide the fit is acceptable |
| Safe to self-serve? | **Yes** | **No** |

**Recommendation: build Job A now as fully self-service, and keep Job B as a
once-a-year controlled refresh with a documented script.** Trying to make Job B a
button is where this gets expensive and risky, and it is not actually what most
users need.

---

## Where each input comes from

| Input | Source | Machine-retrievable? |
|---|---|---|
| CE-QUAL-W2 scenario temperatures | Reclamation temperature modeling team, emailed `.xlsx` | **No.** A human deliverable. Upload is the only path. |
| Observed water temperature | USGS NWIS, stations 11446500 (Hazel) / 11446980 (Watt) | **Yes** — a real public API, already used via `dataRetrieval` |
| Carcass survey detail | SacPAS (Columbia Basin Research, UW) | **Yes** — see below |
| GrandTab adult escapement | SacPAS | **Yes** |
| Water-year hydrologic index (HCI) | SacPAS | **Yes** |
| CWT releases / SAR | RMIS, currently a hand-built spreadsheet | Partly — would need work |
| PHABSIM WUA → carrying capacity | Reclamation | No, and it rarely changes |
| Hydropower revenue by scenario | Internal valuation | No, and it rarely changes |

### The useful discovery about SacPAS

The three SacPAS datasets we depend on are **officially scriptable**. Their query
pages carry this option:

> *"Generate Query Result Link Only — Check box and click 'Submit Query' button to
> generate data link for querying results directly from scripts and automated
> processes."*

So this is **not scraping**. SacPAS intends for automated querying and provides
the mechanism. Someone generates the link once from the web form, we store it, and
the app or a refresh script can pull current data from it.

Sources: [carcass survey detail](https://www.cbr.washington.edu/sacramento/data/query_carcass_detail.html),
[GrandTab escapement](https://www.cbr.washington.edu/sacramento/data/query_adult_grandtab.html),
[SacPAS data home](https://www.cbr.washington.edu/sacramento/data/)

**One caveat that should shape the decision.** Every SacPAS export we hold carries
this notice: *"data presented here are preliminary and subject to revision."* If
the app fetches live at runtime, the model's inputs can change underneath a user
between one session and the next, and two people running "the same" analysis a
month apart get different answers with no record of why. For a decision-support
tool that feeds a published analysis, that is a real problem. It argues for
**pinned snapshots with a deliberate refresh**, not live fetching — regardless of
which option below is chosen.

---

## What we established this week

Two measurements that make Job A much cheaper than we assumed:

1. **The forecast temperature series is exactly periodic.** From 2026 onward every
   projection year has an identical daily profile (verified: 0.000 °C difference
   between years 2030, 2031 and 2075). So egg-to-fry survival only has to be
   computed for **one season per scenario, not 114**.

2. **The expensive part of `precompute.R` is avoidable for this use.** It is slow
   because it loops over 114 years × 36 alternatives × ~400 individual redds. The
   same arithmetic done with cumulative-sum lookups is roughly 100× cheaper.
   Measured: 6 alternatives × 114 years ran in about 14 seconds *including* startup
   and loading a 4.3 MB file. Nine scenarios × one season, in an already-running
   app, is well under a second.

3. **Calibration does not depend on the scenario temperatures.** SAR and
   `rear_surv` are fitted against observed escapement using only the reference
   alternative. New scenario temperatures change the forecast and nothing else.
   *This is what makes Job A safe to self-serve.* It is also exactly why Job B is
   different: new escapement data **does** change the calibration.

---

## The options

### Option 1 — Upload a scenario, explore it, download results *(session-only)*

The user drags in the temperature deliverable exactly as the modeling team sends
it. The app validates it, recomputes the forecast, and shows the new scenarios
beside the published baseline. Nothing is saved server-side; closing the tab
discards it.

- **Effort:** ~2–3 days
- **Runs on:** the existing shinyapps.io deployment, no changes
- **Who can use it:** anyone with the URL, no login, no R
- **Data foundation:** frozen at the current snapshot
- **Risk:** very low — cannot affect anyone else's results or the published numbers

**Good for:** "what would a 2026 met year look like under these nine scenarios?"
**Does not solve:** incorporating new carcass or escapement years.

### Option 2 — Option 1, plus pinned SacPAS snapshots with a visible "as of" date

Same as Option 1, but carcass/escapement/HCI are bundled as dated snapshots pulled
from SacPAS, and every screen and export shows which snapshot is in use.

- **Effort:** +1 day over Option 1
- **Benefit:** removes ambiguity about which vintage of data produced a result;
  makes the app's outputs reproducible and citable
- **Still does not** re-fit anything — the snapshot is used exactly as the current
  data is used

**This is the smallest change that makes the app trustworthy over time,** and I would
do it regardless of what else is chosen.

### Option 3 — Annual refresh script *(recommended companion to 1 + 2)*

A single documented script — `refresh_data_year.R` — that pulls current SacPAS and
NWIS data, re-fits the spawn-timing model and the calibration, regenerates
`app_data`, and writes a short report showing what moved and by how much.

- **Effort:** ~3–4 days, most of it in the comparison report and the checks
- **Run by:** one person, once a year, deliberately
- **Why a script and not a button:** re-calibration is a scientific judgment.
  Someone has to look at the new fit and decide it is acceptable before it becomes
  the official basis for a recommendation. Automating that decision is the actual
  risk here, not the computation.
- **Succession benefit:** this is the thing that removes the single-point-of-failure
  problem. The knowledge stops being in one person's head and becomes a script with
  a checklist.

### Option 4 — Full self-service including new carcass years

Users upload or trigger a fetch of new carcass/escapement data, the app re-fits and
re-calibrates in-session.

- **Effort:** 2–3 weeks
- **Requires:** authentication, an admin role, cloud storage (shinyapps.io wipes
  local disk on restart), version history, a rollback path
- **My honest view: do not do this.** It puts model re-calibration in the hands of
  users who cannot evaluate whether the new fit is sound, and it makes the app's
  official numbers a moving target. The demand for it is also probably low —
  most people want to explore scenarios, not re-fit a life-cycle model.

---

## Recommendation

**Options 1 + 2 + 3.** Roughly a week of work total.

- Anyone can explore a new temperature scenario, unaided, from a URL — which is the
  actual ask.
- Every result is stamped with the data vintage that produced it.
- The annual data refresh becomes a documented, repeatable procedure rather than
  tribal knowledge.
- Model re-calibration stays where it belongs: with a person who can judge the fit.

## Open questions for Brian

1. **Does a new met year *add* a fifth year (45 alternatives, weights 0.2 each) or
   *replace* one for comparison?** This changes the MCDA weighting and the
   consequence table. Currently four met years at 0.25 each.
2. **Who owns the annual refresh** once it is a script, and who signs off that a new
   calibration is acceptable?
3. **Do we want new scenario runs to be citable** — i.e. should the app stamp
   exports with a run ID and data vintage so a number in a memo can be traced back?
4. **Hydropower revenue and steelhead metrics for a new year** — the steelhead
   metric can be computed from the uploaded temperatures, but revenue figures come
   from a separate valuation. Should the app default to the published per-scenario
   values and let the user override?
5. **Is shinyapps.io the long-term home?** The free/basic tier idles out and has
   memory limits. If this becomes a standing tool for the region, a Posit Connect
   or internal server deployment would be steadier.

---

## Appendix — what would actually break if we did nothing

Worth stating plainly, since it is the reason this is on the agenda:

- The temperature deliverable arrives as a spreadsheet that only maps into the model
  through a script one person maintains.
- Adding a year of carcass or escapement data requires re-running a pipeline whose
  slow path takes a long time and whose parameters are set in code.
- `app_data/` is generated output with no recorded provenance — nothing in the files
  says which data vintage or which commit produced them.
- The result is that every request routes through one person, and the institutional
  knowledge is not written down anywhere a successor could pick up.

Options 2 and 3 are aimed squarely at that last point, and they matter more than the
upload feature does.
