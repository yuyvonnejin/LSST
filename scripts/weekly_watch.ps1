# Weekly lsst_watch run: survey + Claude delta report.
# Registered in Windows Task Scheduler (see README, "Scheduled runs").

$repo = "C:\Users\Yvonne Jin\Documents\code\LSST"
$stamp = Get-Date -Format "yyyy-MM-dd"
Set-Location $repo

New-Item -ItemType Directory -Force "$repo\logs" | Out-Null
New-Item -ItemType Directory -Force "$repo\output\history" | Out-Null

# 1. Archive the previous summary (dated by its own mtime)
$summary = "$repo\output\lsst\summary.json"
if (Test-Path $summary) {
    $prevDate = (Get-Item $summary).LastWriteTime.ToString("yyyy-MM-dd")
    Copy-Item $summary "$repo\output\history\lsst_summary_$prevDate.json" -Force
}

# 2. Full-catalog survey run (~2 h)
& "$repo\venv\Scripts\python.exe" -m src.watch --survey lsst *>> "$repo\logs\watch_$stamp.log"

# 3. Headless Claude delta report (read-only; no permissions needed)
$prompt = "You are in the lsst_watch repo. Read output/lsst/summary.json " +
    "and the most recent previous summary in output/history/. Report " +
    "concisely: coverage change (n_ok vs previous), new flare candidates " +
    "(for each new source_id in targets_with_candidates, open its " +
    "output/lsst/target_<id>.json and give distance_pc, band, " +
    "peak_amplitude_mag, peak_flux_ratio, single_epoch, event mjd), and " +
    "whether n_error suggests a broker problem. Plain text, no emoji."
& claude -p $prompt | Out-File "$repo\logs\report_$stamp.md" -Encoding utf8
