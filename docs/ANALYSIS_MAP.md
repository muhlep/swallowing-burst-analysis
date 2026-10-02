# Analysis map

| Manuscript analysis | MATLAB function |
|---|---|
| Common-filter LCMV source reconstruction | `burst_source_reconstruction` |
| Joint sign-deterministic ROI-PC1 extraction | `burst_extract_roi_timeseries` |
| Pooled robust normalization and nested event detection | `burst_detect_events` |
| Mirror-padded analytic amplitude and event reconstruction | `burst_mirror_padded_events` |
| Theta beta count and morphology inference | `burst_event_features` |
| Theta beta and gamma cortical-EMG specificity families | `burst_cortical_emg_overlap` |
| Cropped-envelope 630-edge comparison family | `burst_edge_family_cropped` |
| Primary mirror-padded 630-edge family | `burst_edge_family_mirror_padded` |
| Right-boundary count sensitivity | `burst_boundary_count_sensitivity` |
| Participant-specific forward-inverse spatial leakage | `burst_spatial_leakage_control` |
| Shared-ROI morphology sensitivity | `burst_shared_roi_morphology` |
| Unified marker-relative timing family | `burst_marker_relative_timing` |
| Mirror-padded bidirectional pairwise orthogonalization | `burst_orthogonalized_edges` |
| Bounded onset-lag sensitivity | `burst_bounded_onset_lag` |
| Low-beta and high-beta prevalence inference | `burst_beta_subband_prevalence` |
| Equal-duration beta-subband sensitivity | `burst_beta_subband_duration_matched` |

For the cortical-EMG analysis, call `burst_cortical_emg_overlap` once with
`["theta";"alpha";"beta"]` and once with `["lowgamma";"highgamma"]`.
This retains the two separate multiplicity families described in the methods.

Low-beta and high-beta source time series are reconstructed with
`burst_source_reconstruction` using 13-20 Hz and 20-30 Hz, respectively,
and events are detected with `burst_detect_events`. A separate full-cohort
runner is unnecessary for understanding or reconstructing those steps.

