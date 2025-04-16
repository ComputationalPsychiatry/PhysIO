---
name: PhysIO Bug report
about: Report a bug in PhysIO
title: "[BUG] "
labels: ''
assignees: mrikasper

---

**Describe the bug**
1. Please provide a clear and concise description of what the bug is. 
2. Ideally, provide screenshots of the unexpected behavior, for example
    a. the figure entitled "Preproc: Diagnostics for raw physiological time series" (including the cardiac 
        cycle length curve, respiration histogram), zoomed into time window of concern.
    b. the raw physiological trace data figure entitled "Read-In: Raw Physiological Logfile Data"

**To Reproduce**
Ideally, please provide the following files:
1. Error message in Matlab command window (`.txt` or paste as text in issue)
2. Matlab script or SPM batch file (saved as `.m`) that you ran to call PhysIO
3. Example logfile data corresponding to the Matlab script (`.log`, `.puls`, .`resp` , `.dcm` etc. files), to enable to re-run your script

**Expected behavior**
A clear and concise description of what you expected to happen.

**Software Environment:**
 - OS: [e.g. Windows, MacOS, Linux]
 - PhysIO Version: [e.g. v9.0.1, check via calling `tapas_physio_version()` in the Matlab command window]
- SPM version: none (Matlab standalone) / spm12 / spm25
