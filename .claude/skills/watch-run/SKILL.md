---
name: watch-run
description: Run the lsst_watch alert survey (Rubin/LSST and optionally ZTF), archive the summary, and report what changed since the previous run. Use when the user asks to run the watch, check for new flares, or update the survey.
---

# watch-run: periodic alert-survey run

Run the flare watch over the M dwarf target catalog and report
changes since the last run. Arguments the user may pass: a survey
("lsst", "ztf", or "both"; default lsst) and a target limit
(default: full catalog).

## Steps

1. Archive the previous summary so the new run can be compared:
   for each survey about to run, if `output/<survey>/summary.json`
   exists, copy it to
   `output/history/<survey>_summary_<YYYY-MM-DD>.json`
   (date = file's last-modified date, not today).

2. Refresh the target catalog only if the user asks, or if
   `data/targets.json` is missing:
   `./venv/Scripts/python.exe -m src.targets`

3. Run the watch (full catalog takes ~2 h at 2000 targets; use
   run_in_background for anything over ~100 targets):
   `./venv/Scripts/python.exe -m src.watch --survey lsst`
   Add `--limit N` if the user asked for a bounded run.
   For "both", run lsst first, then ztf.

4. When finished, read `output/<survey>/summary.json` and the most
   recent archived summary for that survey. Report:
   - n_ok / n_no_match / n_error now vs previous (coverage growth)
   - n_flare_candidates now vs previous
   - any new source_ids in targets_with_candidates -- for each new
     one, open its `output/<survey>/target_<id>.json` and show
     distance_pc, band, peak_amplitude_mag, peak_flux_ratio,
     single_epoch, and event mjd
   - errors if n_error rose noticeably (broker outage?)

5. If there are new flare candidates, remind the user these are
   candidates (see README "Reading the output"), and that
   single_epoch events cannot be confirmed from this data alone.

6. Commit the new summary (never the per-target cards, output/ is
   gitignored -- commit nothing unless the user asks; just report).

## Interpreting trends

- n_ok rising run over run = Rubin coverage reaching more of the
  sample. This is expected to climb for months after 2026-06-30.
- A target flipping no_match -> ok means Rubin saw it vary for the
  first time; worth a glance even without a flare candidate.
- If total n_error is high across the run, check the ALeRCE status
  before debugging local code.
