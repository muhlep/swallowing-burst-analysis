# Input schema

## Source reconstruction

The source-reconstruction function receives four study tables:

- `QCLocked`: trial keys and trial-level artifact decisions;
- `LockedSubjectQC`: opaque participant keys and final inclusion decisions;
- `EMGChannelMapping`: dataset locations and EMG channel labels;
- `TimeWindowDefinition`: condition labels and target-window boundaries.

The local subject resolver supplies the head-model, warped-grid, and ROI files.
The completed analysis used 21 predefined ROIs.

## ROI time series

The ROI object contains two conditions, regional arrays with dimensions ROI by
trial by sample, ROI labels, opaque trial keys, time axes, and jointly learned
ROI-PCA weights and centres. Both conditions use the same spatial filters and
ROI weights. Each target window contains 600 samples at 600 Hz.

## Event and inventory tables

Common fields are:

- `Subject`: opaque participant pairing key;
- `Condition`: swallowing `EMG2` or reference `EMG3`;
- `Band` and `Signal`;
- `OriginalTrial` or `TrialIndex`;
- `ZThreshold`, `OnsetSample`, `OffsetSample`, and `DurationSamples`;
- `PeakZ` and suprathreshold area;
- event eligibility and boundary-censoring indicators.

Every function validates the fields it requires. The code retains the
study-specific checks for 73 participants, 21 ROIs, and the prespecified
analysis and correction families.

## Controlled inputs not distributed

Participant-specific MEG and MRI data, spatial filters, forward topographies,
identifier mappings, and local path configuration are not included. Numerical
results are supplied separately through the article Source Data and
Supplementary Tables.

