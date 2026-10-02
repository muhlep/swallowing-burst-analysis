# Burst analysis computational methods

This package contains the MATLAB functions needed to understand and reproduce
the computational analyses reported in the manuscript using appropriately
structured MEG and EMG data. It is a transparent methods archive, not a
one-click application.

The 16 functions cover:

- common-filter source reconstruction and ROI time-series extraction;
- nested event detection and regional event features;
- cortical-EMG and cortico-cortical event-overlap inference;
- mirror-padding, boundary, timing, spatial-leakage, and orthogonalization
  sensitivity analyses;
- the bounded onset-lag and beta-subband analyses.

Raw data, participant mappings, local data-access functions, visualization
code, generated results, and internal workflow files are not included.

## Software

- MATLAB R2025a 25.1.0.2943329
- FieldTrip 20231220
- Signal Processing Toolbox

Each function records or exposes the analysis-specific thresholds, frequency
bands, minimum durations, permutation counts, bootstrap counts, random seeds,
and correction-family definitions used in the study.

## Data interface

The source-reconstruction step expects study tables containing trial-level
quality control, final participant inclusion, dataset and EMG-channel mapping,
and analysis-window definitions. A local subject resolver named `pmd15_cfg`
was used to locate the participant-specific head model, warped source grid, and
ROI definitions. Researchers using another dataset should replace that resolver
while retaining the documented source-analysis settings.

Downstream functions operate on opaque participant keys and the event,
inventory, and ROI time-series structures described in
`docs/INPUT_SCHEMA.md`. No participant data are included here.

## Navigation

- `docs/ANALYSIS_MAP.md` maps manuscript analyses to functions.
- `docs/INPUT_SCHEMA.md` summarizes the expected inputs.
- `docs/SELECTION_RATIONALE.md` records what was deliberately excluded.

Some controlled MAT files use historical field names such as `ROIPilot`.
Those field names are retained only to preserve compatibility with the actual
analysis inputs; the public function names and documentation use descriptive
terminology throughout.

