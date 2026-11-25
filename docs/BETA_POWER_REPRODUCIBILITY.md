# Exploratory beta-band power comparison

`matlab/burst_beta_power.m` computes the paired swallowing versus post-swallow
beta-power comparison on the source-resolved time series and valid trial
support used for cortical event detection.

```matlab
[summary, qc] = burst_beta_power(roiFiles, corticalInventory);
```

- `roiFiles` is a vector of paths to the 73 beta-band ROI MAT files.
  Each contains a `ROIPilot` structure with `subject`, `band`, `bandLabel`,
  `conditions`, `roiLabels`, `roiData`, `time` and `trialInventory`.
  The two data arrays are indexed by ROI, retained trial and time sample.
  The conditions are EMG2 (swallowing) and EMG3 (post-swallow reference).
  The retained-trial table supplies `Condition` and `OriginalTrial` keys.
- `corticalInventory` is the downstream event-detection trial inventory,
  with `Subject`, `Band`, `Condition`, `Signal`, `OriginalTrial`,
  `ZThreshold` and `TrialValid`. Participant keys are opaque pairing
  labels. The function verifies exact ROI-trial key agreement.

For each condition, mean squared amplitude over the 600 target samples is
averaged over valid trials for each ROI. The participant statistic is the
median ROI log ratio (swallowing divided by reference). The mean of those
73 values is back-transformed to percent change. Inference uses a
20,000-resample participant bootstrap for the percentile interval and
100,000 two-sided participant-level sign flips with plus-one correction.
The random seed is 31003.

The reproduction inputs are controlled-access data. This addition includes
the computational method and input description only.
