# Deploying SalmonCountR to shinyapps.io

**The short version.** In a fresh R session, from the repo root:

```r
source("SalmonCountR/deploy.R")
arg_deploy(account = "reclamation-bdo-science")
```

That is the whole thing. No files need to be moved, renamed or rearranged
first. `arg_deploy()` checks the bundle before uploading and refuses to deploy
if the app does not build.

---

## One-time setup

You need `rsconnect` and a configured account token.

```r
install.packages("rsconnect")
```

Then get the token from <https://www.shinyapps.io> → your name (top right) →
**Tokens** → **Show** → **Copy to clipboard**, and paste the
`rsconnect::setAccountInfo(...)` line it gives you into R. Check it took:

```r
rsconnect::accounts()
```

The app has previously been deployed to two accounts, recorded in
`SalmonCountR/rsconnect/shinyapps.io/`:

| Account | URL |
|---|---|
| `reclamation-bdo-science` | the shared Reclamation account |
| `e7ywmx-alex-vaisvil` | Alex's personal account (historical) |

Pass whichever you have access to as `account=`.

---

## What `arg_deploy()` does

1. **Verifies the bundle.** Copies only the files that would be uploaded into a
   temporary directory, then starts a separate R process there and confirms
   `app.R` evaluates and the Shiny app object constructs. This reproduces the
   server's conditions — in particular, no `.git` above the app directory. A
   missing data file or a path that only resolves on your machine fails here,
   in seconds, instead of after a five-minute upload.
2. **Uploads** with an explicit file list.

To check without deploying:

```r
source("SalmonCountR/deploy.R")
arg_verify_bundle()
arg_bundle_files()   # the exact list that would be uploaded
```

---

## Why this file exists

Three things used to make deploying this app a manual chore.

**1. The path bug — fixed in the code, September 2026.** `global.R` built every
data path with `here::here("SalmonCountR", ...)`. `here()` deliberately ignores
the working directory and walks *up* looking for a `.git` / `.Rproj` marker:

| | `here()` resolved to | `here("SalmonCountR", "app_data")` |
|---|---|---|
| Local | repo root (finds `.git`) | ✅ `repo/SalmonCountR/app_data` |
| Deployed | app dir (no marker → falls back to the working directory) | ❌ `appdir/SalmonCountR/app_data` |

Shiny already sets the working directory to the app directory in both cases.
`here()` was the one thing overriding that, so the `SalmonCountR` prefix pointed
one level too deep on the server — hence the reshuffle before every deploy.

Paths now resolve against the directory that actually contains `app_data/`.
See `arg_app_dir()` in `years.R`. **Do not reintroduce `here()` in any file the
app sources** (`app.R`, `global.R`, `years.R`, `functions.R`). It is still
correct in `precompute.R`, which runs from the repo root and is never deployed.

**2. `precompute.R` lives in the app directory.** `rsconnect` scans every `.R`
file in the bundle for `library()` calls, so shipping it made the server install
`furrr`, `ordinal`, `MASS`, `ggridges`, `readxl` and `here` — none of which the
running app uses, and two of which compile from source. It is now excluded.

**3. Bundle size.** 16 of `app_data`'s 25 MB were `precompute.R` intermediates
that no runtime code reads (`sim_future.rds` alone is 13 MB). The deployed
bundle is 23 files, 9.2 MB.

Excluded files are listed in `ARG_DEPLOY_EXCLUDE` in `deploy.R`. It is a
blacklist on purpose: a new data file gets shipped by default, and only the
ones named there are dropped. If the app ever starts reading one of them,
delete it from the list — and note that `arg_verify_bundle()` catches the
mistake either way, because it builds the app from the deploy set alone.

`.rscignore` is generated from that same list by `arg_write_rscignore()`. It
only matters if someone deploys with RStudio's **Publish** button instead of
`arg_deploy()`; re-run that function after editing the list.

---

## Refreshing the data

`precompute.R` writes `app_data/`. Run it from the **repo root**, not from
inside `SalmonCountR/`:

```r
source("SalmonCountR/precompute.R")
```

Then redeploy. `arg_verify_bundle()` will tell you if anything the app needs
failed to write.

## Adding an analysis year

See the header of `years.R`. Short version: drop that year's precompute output
in `app_data/<year>/` and add an entry to `ARG_YEARS`. A year whose files are
not present shows in the selector as "(data not loaded)" rather than erroring,
so a half-finished year cannot break the deployed app.
