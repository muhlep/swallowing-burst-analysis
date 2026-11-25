function [Summary, QC] = burst_beta_power(roiFiles, inventory)
% Paired beta power on the trial support used for event detection.
% roiFiles holds paths to the study-format ROI MAT files, in any order.
% inventory is the event detector's cortical trial inventory table.
paths = sort(string(roiFiles(:)));
assert(numel(paths)==73 && all(isfile(paths)), 'Expected 73 existing beta ROI MAT files.');

keep = strcmpi(string(inventory.Band), "beta") & ...
    double(inventory.ZThreshold) == 2 & logical(inventory.TrialValid);
inventory = inventory(keep, {'Subject','Condition','Signal','OriginalTrial'});
lockedSubjects = sort(unique(string(inventory.Subject)));
nSubjects = numel(lockedSubjects);
assert(nSubjects == 73, 'Expected 73 participants in final beta inventory.');
conditions = ["EMG2"; "EMG3"];
assert(isequal(sort(unique(string(inventory.Condition))), conditions), ...
    'Unexpected conditions in final beta inventory.');

logRatio = nan(nSubjects, 1);
trialCounts = zeros(nSubjects, 2);
seen = false(nSubjects, 1);
for fileIndex = 1:numel(paths)
    loaded = load(char(paths(fileIndex)), 'ROIPilot');
    roi = loaded.ROIPilot;
    assert(strcmpi(string(roi.bandLabel), "beta"), 'Unexpected band label.');
    assert(isequal(double(roi.band(:))', [13 30]), 'Expected 13-30 Hz.');
    subjectIndex = find(lockedSubjects == string(roi.subject));
    assert(isscalar(subjectIndex), 'ROI participant absent from final event-analysis.');
    assert(~seen(subjectIndex), 'Duplicate ROI participant at file %d.', fileIndex);
    seen(subjectIndex) = true;

    labels = string(roi.roiLabels(:));
    assert(numel(labels) == 21 && numel(unique(labels)) == 21, ...
        'Expected 21 unique ROI labels.');
    assert(isequal(string(roi.conditions(:)), conditions), ...
        'Unexpected ROI condition order.');
    participantInventory = inventory(string(inventory.Subject) == ...
        lockedSubjects(subjectIndex), :);
    powerByROI = nan(21, 2);
    for c = 1:2
        x = double(roi.roiData{c});
        nTrials = size(x, 2);
        assert(size(x, 1) == 21 && nTrials > 0 && size(x, 3) == 600 && ...
            all(isfinite(x), 'all'), 'Invalid ROI array at file %d.', fileIndex);
        savedTime = double(roi.time{c}(:));
        assert(numel(savedTime) == 600 && ...
            all(abs(diff(savedTime) - 1/600) < 1e-8), ...
            'Expected 600 samples at 600 Hz.');

        sourceInventory = roi.trialInventory;
        sourceTrials = double(sourceInventory.OriginalTrial( ...
            string(sourceInventory.Condition) == conditions(c)));
        sourceTrials = sourceTrials(:);
        assert(numel(sourceTrials) == nTrials && all(isfinite(sourceTrials)) && ...
            numel(unique(sourceTrials)) == nTrials, 'Invalid saved trial keys.');

        % Check exact condition, ROI and original-trial keys against detection.
        q = string(participantInventory.Condition) == conditions(c);
        actual = string(participantInventory.Signal(q)) + "|" + ...
            string(participantInventory.OriginalTrial(q));
        expected = repelem(labels, nTrials, 1) + "|" + ...
            repmat(string(sourceTrials), numel(labels), 1);
        assert(isequal(sort(actual(:)), sort(expected(:))), ...
            'Valid ROI-trial keys differ from the event inventory at file %d, condition %s.', ...
            fileIndex, char(conditions(c)));
        trialCounts(subjectIndex, c) = nTrials;

        % Average squared samples within each trial, then average valid trials.
        powerByROI(:, c) = mean(mean(x.^2, 3), 2);
    end
    assert(all(isfinite(powerByROI) & powerByROI > 0, 'all'), ...
        'Non-positive or non-finite power at file %d.', fileIndex);
    logRatio(subjectIndex) = median(log(powerByROI(:, 1)) - log(powerByROI(:, 2)));
end
assert(all(seen) && all(isfinite(logRatio)), 'Incomplete participant pairing.');

seed = 31003;
nBootstrap = 20000;
nPermutations = 100000;
originalRng = rng;
restoreRng = onCleanup(@() rng(originalRng)); %#ok<NASGU>
rng(seed, 'twister');
observed = mean(logRatio);
bootstrapMeans = zeros(nBootstrap, 1);
moreExtreme = 0;
batchSize = 2000;
for startIndex = 1:batchSize:nBootstrap
    stopIndex = min(startIndex + batchSize - 1, nBootstrap);
    index = randi(nSubjects, nSubjects, stopIndex - startIndex + 1);
    bootstrapMeans(startIndex:stopIndex) = mean(logRatio(index), 1)';
end
ci = interp1((1:nBootstrap)', sort(bootstrapMeans), ...
    1 + (nBootstrap - 1)*[0.025 0.975], 'linear');
for startIndex = 1:batchSize:nPermutations
    nThisBatch = min(batchSize, nPermutations - startIndex + 1);
    signs = 2*(rand(nSubjects, nThisBatch) >= 0.5) - 1;
    nullMeans = (logRatio' * signs) / nSubjects;
    moreExtreme = moreExtreme + sum(abs(nullMeans) >= abs(observed) - 8*eps(abs(observed)));
end
pValue = (1 + moreExtreme)/(1 + nPermutations);
percent = @(x) 100*expm1(x);
Summary = table(nSubjects, observed, ci(1), ci(2), exp(observed), ...
    percent(observed), percent(ci(1)), percent(ci(2)), pValue, ...
    nBootstrap, nPermutations, seed, ...
    'VariableNames', {'NPaired','MeanLogPowerRatio','LogCI95Lower','LogCI95Upper', ...
    'GeometricMeanPowerRatio','PowerChangePercent','PercentCI95Lower', ...
    'PercentCI95Upper','TwoSidedSignFlipP','NBootstrap','NPermutations','Seed'});
QC = table(nSubjects, 21, 600, 600, sum(trialCounts(:, 1)), ...
    sum(trialCounts(:, 2)), height(inventory), all(seen), ...
    'VariableNames', {'NParticipants','NROIs','SamplesPerTrial','SamplingRateHz', ...
    'NEMG2Trials','NEMG3Trials','NMatchedValidROITrials','ExactTrialSupportMatch'});

end
