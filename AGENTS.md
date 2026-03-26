# AGENTS.md — AI Agent Guide for the PhysIO Toolbox

This file provides context and instructions for AI coding agents (Claude Code,
GitHub Copilot, Cursor, etc.) working in this repository. Read it before making
any changes.

---

## What is PhysIO?

PhysIO is a **MATLAB toolbox** for model-based physiological noise correction of
fMRI data. It takes peripheral physiological recordings (cardiac pulsation via
ECG or pulse oximetry; respiration via breathing belt) as input and produces
nuisance regressor text files as output. These regressors can be fed into any
major fMRI analysis package (SPM, FSL, AFNI, BrainVoyager, etc.) as columns of
a General Linear Model (GLM).

**Do not assume this is a Python project.** All core code is MATLAB (`.m` files).
There is no `setup.py`, no `pyproject.toml`, no `package.json`.

**Part of TAPAS.** PhysIO is the physiological noise correction component of the
[TAPAS software collection](https://www.translationalneuromodeling.org/tapas)
(Translational Algorithms for Psychiatry-Advancing Science) by the
Translational Neuromodeling Unit (TNU), University of Zurich and ETH Zurich, now hosted
under the [ComputationalPsychiatry GitHub organization](https://github.com/ComputationalPsychiatry).

---

## Repository Structure

```
PhysIO/
├── code/                     # All MATLAB source code
│   ├── tapas_physio_new.m    # *** MASTER PARAMETER DEFINITION FILE ***
│   │                         # Defines the `physio` struct with inline comments
│   │                         # on every parameter. Read this first to understand
│   │                         # the data model.
│   ├── tapas_physio_main_create_regressors.m  # Top-level pipeline entry point
│   ├── tapas_physio_init.m   # Sets up MATLAB path and SPM toolbox link
│   ├── tapas_physio_review.m # Recreates QA figures from saved physio.mat
│   ├── tapas_physio_report_contrasts.m  # SPM F-test/tSNR QA reporting
│   ├── readin/               # Format-specific log-file readers
│   │   ├── tapas_physio_read_physlogfiles_bids.m
│   │   ├── tapas_physio_read_physlogfiles_biopac_mat.m
│   │   ├── tapas_physio_read_physlogfiles_ge.m
│   │   ├── tapas_physio_read_physlogfiles_philips.m
│   │   ├── tapas_physio_read_physlogfiles_siemens.m  (VB format, ideaCmdTool)
│   │   ├── tapas_physio_read_physlogfiles_siemens_tics.m  (VD/VE/Tics format, XA60, CMRR)
│   │   └── tapas_physio_read_physlogfiles_custom.m
│   ├── preproc/              # Preprocessing: peak detection, filtering
│   ├── model/                # Noise model creation (RETROICOR, HRV, RVT, etc.)
│   ├── assess/               # Model performance assessment
│   └── utils/                # Utility functions
├── examples/                 # Example scripts and input logfiles from different vendors, not part of this repository (see note below)
│   ├── BIDS/
│   ├── GE/
│   ├── Philips/
│   ├── Siemens_VB/
│   ├── Siemens_VD/
│   └── HCP/
├── tests/
│   ├── unit/                 # Unit tests (MATLAB testing framework)
│   └── integration/          # Integration tests for all examples
└── test-reference-results/   # Golden reference output files, not part of this repository (see note below)
├── docs/                     # Additional documentation
├── README.md
├── CHANGELOG.md
└── AGENTS.md                 # This file
```

> **Note:** Example data (logfiles) in examples/ subfolder are **not** stored in this repository. They
> are downloaded separately by calling `tapas_physio_download_example_data()` from
> the MATLAB command line after installation or can be cloned directly from
>  https://github.com/ComputationalPsychiatry/PhysIO-Examples

> **Note:**  Test reference results in test-reference-resuts/ subfolder are required for running unit 
> and integration tests in the tests/ subfolder, but **not** stored in this repository. They
> are downloaded separately by calling `tapas_physio_download_test_reference_results()` from
> the MATLAB command line after installation or can be cloned directly from 
> https://github.com/ComputationalPsychiatry/PhysIO-Test-Reference-Results

---

## The Central Data Structure: `physio`

Everything in PhysIO revolves around a single MATLAB struct called `physio`.
It is constructed by `tapas_physio_new()` and populated/passed through the whole
pipeline.

**Before modifying any pipeline code, open and read `code/tapas_physio_new.m`.**
Every parameter is defined there with inline comments explaining its purpose and
valid values. This is the authoritative technical reference — it supersedes the
scientific papers, which may use older variable names.

Key top-level substructures of `physio`:

| Substructure | Purpose |
|---|---|
| `physio.log_files` | Input file paths and vendor/format selection |
| `physio.scan_timing` | fMRI acquisition parameters (TR, nSlices, nVols, etc.) |
| `physio.preproc` | Preprocessing options (cardiac/respiratory filtering, peak detection) |
| `physio.model` | Noise model selection and order (RETROICOR, HRV, RVT, noise ROIs, motion) |
| `physio.model.retroicor` | RETROICOR Fourier expansion orders |
| `physio.model.hrv` | Heart rate variability model and delays |
| `physio.model.rvt` | Respiratory volume per time model and delays |
| `physio.model.noise_rois` | aCompCor-style noise ROI regressors (requires SPM) |
| `physio.model.other` | Include additional regressor text file |
| `physio.model.movement` | Motion parameter regressors (6/12/24 parameters) |
| `physio.model.R` | **Output:** Combined confound regressor matrix from all selected models, to be put into GLM for 1st level statistical parametric mapping analysis in any fMRI analysis package |
| `physio.ons_secs` | **Output:** intermediate time courses (set during pipeline execution) |
| `physio.verbose` | Logging verbosity and figure saving settings |
| `physio.save_dir` | Output directory |

The `physio.ons_secs` substructure holds all intermediate results and is saved
to `physio.mat` in the output directory after a run. Key fields include:
- `ons_secs.t` — time vector
- `ons_secs.c` / `ons_secs.r` — raw cardiac / respiratory waveforms
- `ons_secs.cpulse` — detected cardiac pulse event times (R-peaks)
- `ons_secs.hr` / `ons_secs.rvt` — heart rate / respiratory volume per time (one per scan volume)
- `ons_secs.fr` — filtered respiratory trace

---

## Function Naming Conventions

All toolbox functions are prefixed `tapas_physio_`. Do not create functions
without this prefix. Examples:

- `tapas_physio_new` — constructor
- `tapas_physio_main_create_regressors` — main pipeline
- `tapas_physio_read_physlogfiles_<vendor>` — vendor-specific readers
- `tapas_physio_preproc_cardiac` / `tapas_physio_preproc_respiratory` — preprocessing
- `tapas_physio_create_retroicor_regressors` — RETROICOR model
- `tapas_physio_create_rvt_regressors` — RVT model
- `tapas_physio_create_hrv_regressors` — HRV model
- `tapas_physio_review` — QA figure recreation

---

## How to Run PhysIO

### Minimal example (MATLAB only, no SPM)

```matlab
% 1. Initialize (adds code/ to path, links SPM toolbox folder if SPM found)
tapas_physio_init();

% 2. Create a default physio struct
physio = tapas_physio_new();

% 3. Set required parameters (see tapas_physio_new.m for all options)
physio.log_files.vendor           = 'Philips';   % or 'BIDS', 'Siemens', 'GE', 'Biopac_Mat', 'Custom', etc.
physio.log_files.cardiac          = 'SCANPHYSLOG.log';
physio.log_files.respiratory      = 'SCANPHYSLOG.log';
physio.log_files.scan_timing      = [];
physio.scan_timing.sqpar.Nslices  = 37;
physio.scan_timing.sqpar.TR       = 2.5;
physio.scan_timing.sqpar.Nscans   = 200;
physio.model.retroicor.order.c    = 3;
physio.model.retroicor.order.r    = 4;
physio.model.retroicor.order.cr   = 1;
physio.save_dir                   = './output';

% 4. Run the pipeline
physio = tapas_physio_main_create_regressors(physio);
```

### Via SPM Batch Editor (with SPM)

1. Run `tapas_physio_init()` from the MATLAB command line.
2. Start SPM: `spm fmri`
3. Open the Batch Editor and navigate to: **SPM → Tools → TAPAS PhysIO Toolbox**
4. Or load an existing batch file (`*_spm_job.m` / `*_spm_job.mat`) from any
   `examples/` subfolder.

### Using example scripts

```matlab
tapas_physio_init();
cd('examples/Philips/ECG3T');
% Run example_philips_ecg3t_matlab_script.m (Matlab-only)
% or example_philips_ecg3t_spm_job.m (SPM Batch Editor)
run('example_philips_ecg3t_matlab_script.m');
```

---

## Supported Input Formats (Vendors)

The `physio.log_files.vendor` parameter selects the reader. Valid string values
and their corresponding formats:

| Vendor string | Format | Key files |
|---|---|---|
| `'BIDS'` | Brain Imaging Data Structure `.tsv[.gz]` + `.json` sidecar | `*_physio.tsv.gz`, `*_physio.json` |
| `'Biopac_Mat'` | BioPac `.mat` export | `.mat` (variables: `data`, `isi`, `isi_units`, `labels`, `start_sample`, `units`) |
| `'Biopac_Txt'` | BioPac `.txt` export | 4-column text (resp, GSR, cardiac, trigger) |
| `'GE'` | General Electric | `ECGData_*`, `RespData_*` (one amplitude per line) |
| `'Philips'` | Philips SCANPHYSLOG | `SCANPHYSLOG<DateTime>.log` |
| `'Siemens'` | Siemens VB (manual/IdeaCmdTool) | `.ecg`, `.resp`, `.puls` + optional `.dcm` |
| `'Siemens_Tics'` | Siemens VD/VE/XA (CMRR) | `*_Info.log`, `*_PULS.log`, `*_RESP.log`, `*_ECG.log` |
| `'Siemens_HCP'` | Human Connectome Project | `*_Physio_log.txt` |
| `'Custom'` | Any vendor (plain text) | Separate cardiac and respiratory text files, one amplitude per line |

The corresponding MATLAB reader for each vendor lives in
`code/readin/tapas_physio_read_physlogfiles_<vendor_lowercase>.m`.

---

## Noise Models

PhysIO implements several additive noise models. All are controlled via
`physio.model.*` and output columns in a single multiple-regressors text file.

**Regressor column order in output file:**

1. RETROICOR cardiac [2 × `order.c`]
2. RETROICOR respiratory [2 × `order.r`]
3. RETROICOR cardiac × respiratory interactions [4 × `order.cr`]
4. HRV [`nDelaysHRV`]
5. RVT [`nDelaysRVT`]
6. Noise ROIs / aCompCor [nROIs × (nComponents + 1)]
7. Other (additional text file columns)
8. Motion [6, 12, or 24 depending on `model.movement.type`]

**RETROICOR** (Glover 2000; Harvey 2008): Fourier expansion of estimated cardiac
and respiratory phases. Standard order: cardiac 3rd, respiratory 4th,
cardio-respiratory interaction 1st (giving 18 regressors total).

**RVT** (Birn 2006; Harrison et al. 2021): Respiratory volume per time with a
response function convolution. For respiratory trace preprocessing and RVT, the
Hilbert-based method (Harrison et al. 2021) is used; please cite that paper when
using RVT regressors or respiratory RETROICOR.

**HRV** (Chang 2009): Heart rate variability with a cardiac response function
convolution.

**Noise ROIs / aCompCor** (Behzadi 2007): PCA of signal from noise regions of
interest (requires NIfTI image read via SPM).

**Motion** (Siegel 2014, Friston 1997): 6, 12 (+ squared), or 24 (+ derivatives + squared)
motion parameter regressors.

---

## Testing

Tests use the **MATLAB unit testing framework** (`matlab.unittest`).

```matlab
% Run all unit tests
cd tests/unit
results = runtests('.');

% Run all integration tests (requires downloaded example data)
cd tests/integration
results = runtests('.');
```

Integration tests compare pipeline output against reference data stored in
`test-reference-results/` subfolder (from different repository, see note above). 
When adding a new example or fixing a bug that legitimately changes output, 
update the reference results accordingly.

Both SPM Batch Editor and MATLAB-only code paths are tested for all examples.

---

## SPM Integration Details

PhysIO works standalone in MATLAB, but gains extra capabilities when SPM12 or SPM25 is
installed:

- **Batch Editor GUI** — full graphical interface under SPM → Tools → TAPAS PhysIO Toolbox
- **Pipeline dependencies** — automatic input of realignment parameters, feed-in
  to SPM's multiple-regressors GLM slot
- **Model assessment** — F-contrast maps and tSNR reports via
  `tapas_physio_report_contrasts()`
- **Noise ROIs** — reading NIfTI mask images via `spm_vol` / `spm_read_vols`

To install for SPM: copy or symlink the `PhysIO/` folder into
`<path-to-spm12>/toolbox/PhysIO`. `tapas_physio_init()` handles this check and
will warn if the link is missing.

**Important path note:** Never use `addpath(genpath(...))` for SPM — this causes
function name conflicts. Always use `addpath('<path-to-spm12>')` (no `genpath`).

---

## No-Matlab Deployment

PhysIO can be run without a MATLAB license via:
- **[NeuroDesk](https://www.neurodesk.org/)** — containerized SPM+PhysIO
- **[CBRAIN](https://mcin.ca/technology/cbrain/)** — web-based processing

These use a compiled MATLAB Runtime version. Developers wishing to compile SPM
with PhysIO can use `spm_make_standalone`; the PhysIO codebase is compatible
with the MATLAB Compiler.

---

## How to Cite

All publications using PhysIO must cite:

1. **Core PhysIO paper:** Kasper et al. (2017). *The PhysIO Toolbox for Modeling
   Physiological Noise in fMRI Data.* Journal of Neuroscience Methods 276, 56–72.
   https://doi.org/10.1016/j.jneumeth.2016.10.019

2. **If using RVT or respiratory RETROICOR preprocessing:** Harrison et al. (2021).
   *A Hilbert-based method for processing respiratory timeseries.*
   NeuroImage, 117787. https://doi.org/10.1016/j.neuroimage.2021.117787

3. **TAPAS collection:** Frässle et al. (2021). *TAPAS: an open-source software
   package for Translational Neuromodeling and Computational Psychiatry.*
   Frontiers in Psychiatry 12, 857.
   https://doi.org/10.3389/fpsyt.2021.680811


Methods section snippet:

> The analysis was performed using the Matlab PhysIO Toolbox ([1,2], version x.y.z,
> open-source code available as part of the TAPAS software collection: [3],
> https://www.translationalneuromodeling.org/tapas)

---

## Contributing

- New contributors must add themselves to the **Contributor License Agreement** as
  part of their first pull request (referenced in README.md), because of the GPLv3 license requirements.
- Follow the existing `tapas_physio_` function prefix convention.
- Add/update integration tests and reference results for any new functionality or
  bug fixes that change pipeline output.
- Support questions go to [GitHub Issues](https://github.com/ComputationalPsychiatry/PhysIO/issues).
  Older discussions (2018–2025) are in the
  [previous TAPAS repository issues](https://github.com/translationalneuromodeling/tapas/issues).

---

## Key Documentation Links

- **Parameter reference:** `edit tapas_physio_new` in MATLAB, or `help tapas_physio_new`
- **Wiki (most up-to-date user guide):** https://github.com/ComputationalPsychiatry/PhysIO/wiki
  - [QUICKSTART](https://github.com/ComputationalPsychiatry/PhysIO/wiki/QUICKSTART)
  - [FAQ](https://github.com/ComputationalPsychiatry/PhysIO/wiki/FAQ)
  - [Supported logfile formats](https://github.com/ComputationalPsychiatry/PhysIO/wiki/MANUAL_PART_READIN)
  - [Example datasets](https://github.com/ComputationalPsychiatry/PhysIO/wiki/EXAMPLES)
- **CHANGELOG:** `CHANGELOG.md`
- **Scientific background:** https://doi.org/10.1016/j.jneumeth.2016.10.019
