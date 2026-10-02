function Results = burst_mirror_padded_events(varargin)
% Bounded sensitivity for Hilbert edge effects in the locked healthy cohort.
%
% This function keeps the locked participants, retained trials, source
% reconstruction, ROI-PC weights, bands, thresholds and duration rules
% unchanged. The only analytical change is calculation of cortical analytic
% amplitude after symmetric reflection padding of each 600-sample regional
% time series. Padding is removed before normalization and event detection.
%
% The run is resumable at participant level and never overwrites inputs.

p = inputParser;
addParameter(p,'sourceRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'legacyBurstFile','',@(x)ischar(x)||isstring(x));
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'resumeRunDir','',@(x)ischar(x)||isstring(x));
addParameter(p,'mirrorPadSamples',599,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun,'Set confirmRun=true only after checking all paths.');

sourceRoot = char(string(p.Results.sourceRoot));
legacyBurstFile = char(string(p.Results.legacyBurstFile));
outputRoot = char(string(p.Results.outputRoot));
assert(isfolder(sourceRoot),'Source root missing: %s',sourceRoot);
assert(isfile(legacyBurstFile),'Legacy burst file missing: %s',legacyBurstFile);
assert(isfolder(outputRoot),'Output root missing: %s',outputRoot);

Legacy = loadBurstObject(legacyBurstFile);
subjects = sort(unique(string(Legacy.CorticalInventory.Subject)));
bands = ["theta";"alpha";"beta"];
bandHz = [4 8;8 13;13 30];
conditions = ["EMG2";"EMG3"];
thresholds = [1.5 2 2.5];
fs = 600;
assert(numel(subjects)==73,'Expected 73 healthy participants.');

allFiles = dir(fullfile(sourceRoot,'**','*_ROI_PCA_PILOT.mat'));
allPaths = strings(numel(allFiles),1);
for k = 1:numel(allFiles)
    allPaths(k) = string(fullfile(allFiles(k).folder,allFiles(k).name));
end
standardMask = false(size(allPaths));
for b = 1:numel(bands)
    standardMask = standardMask | contains(upper(allPaths),"_"+upper(bands(b))+"_");
end
standardPaths = allPaths(standardMask);
assert(numel(standardPaths)==73*3, ...
    'Expected exactly 219 standard-band ROI-PCA files; found %d.',numel(standardPaths));

resumeRunDir = char(string(p.Results.resumeRunDir));
if isempty(resumeRunDir)
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    runID = ['BURST_MIRROR_PADDED_EVENTS_' stamp];
    runDir = fullfile(outputRoot,runID);
    assert(~isfolder(runDir),'Refusing overwrite: %s',runDir);
    [ok,msg] = mkdir(runDir); assert(ok,'Cannot create run directory: %s',msg);
    Status = table(subjects,repmat("PENDING",73,1),repmat("",73,1), ...
        nan(73,1),'VariableNames',{'Subject','State','Message','ElapsedSeconds'});
    cfg = struct('runID',runID,'runDir',string(runDir), ...
        'sourceRoot',string(sourceRoot),'legacyBurstFile',string(legacyBurstFile), ...
        'subjects',subjects,'bands',bands,'bandHz',bandHz,'conditions',conditions, ...
        'thresholds',thresholds,'samplingRateHz',fs, ...
        'mirrorPadSamples',double(p.Results.mirrorPadSamples), ...
        'changeFromLockedAnalysis', ...
        "cortical analytic amplitude calculated after symmetric reflection padding");
else
    runDir = resumeRunDir;
    assert(isfolder(runDir),'Resume directory missing: %s',runDir);
    stateHits = dir(fullfile(runDir,'*_STATE.mat'));
    assert(numel(stateHits)==1,'Expected exactly one saved state file.');
    S = load(fullfile(stateHits(1).folder,stateHits(1).name),'cfg','Status');
    cfg = S.cfg; Status = S.Status; runID = char(cfg.runID);
    assert(string(cfg.sourceRoot)==string(sourceRoot),'Resume sourceRoot differs.');
    assert(string(cfg.legacyBurstFile)==string(legacyBurstFile), ...
        'Resume legacyBurstFile differs.');
end
stateFile = fullfile(runDir,[runID '_STATE.mat']);
save(stateFile,'cfg','Status','-v7.3');

timerAll = tic;
for s = 1:numel(subjects)
    sid = subjects(s);
    subjectFile = fullfile(runDir,sprintf('%s_SUBJECT_%s.mat',runID,char(sid)));
    if isfile(subjectFile)
        Status.State(s) = "COMPLETE";
        Status.Message(s) = "RESUMED_FROM_SUBJECT_FILE";
        fprintf('[%d/73] %s already complete; skipped.\n',s,sid);
        continue
    end
    fprintf('\n[%d/73] MIRROR-PADDED HILBERT: %s\n',s,sid);
    subjectTimer = tic;
    try
        eventParts = cell(3,1); inventoryParts = cell(3,1); normParts = cell(3,1);
        inputRows = cell(3,6);
        for b = 1:3
            band = bands(b);
            candidates = standardPaths(contains(standardPaths,sid,'IgnoreCase',true) & ...
                contains(upper(standardPaths),"_"+upper(band)+"_"));
            assert(numel(candidates)==1, ...
                'Expected one ROI-PCA file for %s %s; found %d.',sid,band,numel(candidates));
            roiFile = char(candidates(1));
            R = load(roiFile,'ROIPilot');
            P = R.ROIPilot;
            assert(string(P.subject)==sid && string(P.bandLabel)==band, ...
                'ROI-PCA identity mismatch.');
            assert(isequal(string(P.conditions(:)),conditions), ...
                'Condition order differs from EMG2/EMG3.');
            assert(numel(P.roiLabels)==21,'Expected 21 ROIs.');
            assert(size(P.roiData{1},3)==600 && size(P.roiData{2},3)==600, ...
                'Expected 600 target samples.');
            trialKeys = getTrialKeys(P,conditions);
            minimumSamples = ceil(fs*2/mean(bandHz(b,:)));
            [E,I,N] = detectMirrorPaddedEvents(P.roiData,conditions,trialKeys, ...
                P.roiLabels,band,fs,thresholds,minimumSamples, ...
                double(p.Results.mirrorPadSamples));
            E.Subject = repmat(sid,height(E),1);
            I.Subject = repmat(sid,height(I),1);
            N.Subject = repmat(sid,height(N),1);
            eventParts{b} = E; inventoryParts{b} = I; normParts{b} = N;
            inputRows(b,:) = {sid,band,string(roiFile),size(P.roiData{1},2), ...
                size(P.roiData{2},2),minimumSamples};
            fprintf('  %s | EMG2=%d trials | EMG3=%d trials | events=%d\n', ...
                band,size(P.roiData{1},2),size(P.roiData{2},2),height(E));
            clear R P E I N
        end
        SubjectEvents = vertcat(eventParts{:});
        SubjectInventory = vertcat(inventoryParts{:});
        SubjectNormalization = vertcat(normParts{:});
        SubjectInputs = cell2table(inputRows,'VariableNames', ...
            {'Subject','Band','ROIFile','NTrialsEMG2','NTrialsEMG3', ...
            'MinimumDurationSamples'});
        save(subjectFile,'SubjectEvents','SubjectInventory', ...
            'SubjectNormalization','SubjectInputs','-v7.3');
        Status.State(s) = "COMPLETE"; Status.Message(s) = "OK";
        Status.ElapsedSeconds(s) = toc(subjectTimer);
        save(stateFile,'cfg','Status','-v7.3');
        writetable(Status,fullfile(runDir,[runID '_STATUS.csv']));
        fprintf('%s complete | %.1f min | total %s\n',sid, ...
            Status.ElapsedSeconds(s)/60, ...
            char(duration(0,0,toc(timerAll),'Format','hh:mm:ss')));
    catch ME
        Status.State(s) = "FAILED"; Status.Message(s) = string(ME.message);
        Status.ElapsedSeconds(s) = toc(subjectTimer);
        save(stateFile,'cfg','Status','ME','-v7.3');
        writetable(Status,fullfile(runDir,[runID '_STATUS.csv']));
        rethrow(ME)
    end
end
assert(all(Status.State=="COMPLETE"),'Not all participants completed.');

eventParts = cell(73,1); inventoryParts = cell(73,1);
normParts = cell(73,1); inputParts = cell(73,1);
for s = 1:73
    sid = subjects(s);
    subjectFile = fullfile(runDir,sprintf('%s_SUBJECT_%s.mat',runID,char(sid)));
    S = load(subjectFile,'SubjectEvents','SubjectInventory', ...
        'SubjectNormalization','SubjectInputs');
    eventParts{s} = S.SubjectEvents; inventoryParts{s} = S.SubjectInventory;
    normParts{s} = S.SubjectNormalization; inputParts{s} = S.SubjectInputs;
end
CorticalEvents = vertcat(eventParts{:});
CorticalInventory = vertcat(inventoryParts{:});
CorticalNormalization = vertcat(normParts{:});
InputFiles = vertcat(inputParts{:});
CorticalEvents.GlobalEventID = (1:height(CorticalEvents))';

NestingQC = checkNesting(CorticalInventory);
assert(all(NestingQC.NViolationsInclusive==0) && ...
    all(NestingQC.NViolationsPrimary==0),'Threshold nesting failed.');
LegacyComparison = compareLegacyCounts(Legacy.CorticalEvents, ...
    Legacy.CorticalInventory,CorticalEvents,CorticalInventory);

BurstMirrorPadded = struct( ...
    'CorticalEvents',CorticalEvents, ...
    'CorticalInventory',CorticalInventory, ...
    'CorticalNormalization',CorticalNormalization, ...
    'EMGEvents',Legacy.EMGEvents, ...
    'EMGInventory',Legacy.EMGInventory, ...
    'EMGNormalization',Legacy.EMGNormalization, ...
    'NestingQC',NestingQC, ...
    'LegacyCountComparison',LegacyComparison, ...
    'InputFiles',InputFiles, ...
    'cfg',cfg);

burstOutputFile = fullfile(runDir,[runID '_BURST_MIRROR_PADDED.mat']);
save(burstOutputFile,'BurstMirrorPadded','cfg','-v7.3');
writetable(NestingQC,fullfile(runDir,'NESTING_QC.csv'));
writetable(LegacyComparison,fullfile(runDir,'LEGACY_EVENT_COUNT_COMPARISON.csv'));
writetable(InputFiles,fullfile(runDir,'INPUT_FILE_AUDIT.csv'));
writetable(Status,fullfile(runDir,[runID '_STATUS.csv']));

InputQC = table(73,21,3,600,double(p.Results.mirrorPadSamples), ...
    height(CorticalEvents),height(CorticalInventory),all(Status.State=="COMPLETE"), ...
    'VariableNames',{'NSubjects','NROIs','NBands','TargetSamples', ...
    'MirrorPadSamples','NCorticalEvents','NCorticalInventoryRows','Complete'});
writetable(InputQC,fullfile(runDir,'INPUT_QC.csv'));
Results = struct('InputQC',InputQC,'NestingQC',NestingQC, ...
    'LegacyCountComparison',LegacyComparison,'Status',Status, ...
    'burstOutputFile',string(burstOutputFile),'cfg',cfg);
save(fullfile(runDir,[runID '_RESULTS.mat']),'Results','-v7.3');

fprintf('\nMirror-padded event analysis complete\n%s\n',runDir);
disp(InputQC)
disp(LegacyComparison)
fprintf('BURST OUTPUT\n%s\n',burstOutputFile);
end

function [Events,Inventory,Normalization] = detectMirrorPaddedEvents( ...
    X,conditions,trialKeys,labels,bandLabel,fs,thresholds,minSamples,padRequested)
labels = string(labels(:)); conditions = string(conditions(:));
assert(numel(X)==2 && numel(conditions)==2 && numel(trialKeys)==2, ...
    'Two conditions required.');
assert(size(X{1},1)==numel(labels) && size(X{2},1)==numel(labels), ...
    'Label count mismatch.');
assert(size(X{1},3)==600 && size(X{2},3)==600,'Expected 600 samples.');

Env = cell(2,1);
for c = 1:2
    Env{c} = nan(size(X{c}));
    for r = 1:numel(labels)
        for t = 1:size(X{c},2)
            raw = reshape(double(X{c}(r,t,:)),1,[]);
            assert(all(isfinite(raw)),'Non-finite regional time series.');
            n = numel(raw); padN = min(round(padRequested),n-1);
            padded = [fliplr(raw(1:padN)) raw fliplr(raw(end-padN+1:end))];
            envelopePadded = analyticEnvelopeFFT(padded);
            envelope = envelopePadded(padN+(1:n));
            Env{c}(r,t,:) = reshape(envelope,1,1,n);
        end
    end
end

normRows = cell(numel(labels),8);
eventRows = cell(0,19); inventoryRows = cell(0,16); eventID = 0;
for r = 1:numel(labels)
    pooled = [reshape(double(Env{1}(r,:,:)),[],1); ...
        reshape(double(Env{2}(r,:,:)),[],1)];
    center = median(pooled,'omitnan');
    scale = 1.4826*median(abs(pooled-center),'omitnan');
    fallback = false;
    if ~isfinite(scale) || scale<=eps(max(1,abs(center)))
        scale = std(pooled,0,'omitnan'); fallback = true;
    end
    assert(isfinite(center)&&isfinite(scale)&&scale>0, ...
        'Invalid normalization for %s.',labels(r));
    normRows(r,:) = {labels(r),string(bandLabel),center,scale,numel(pooled), ...
        sum(isfinite(pooled)),fallback,"CORTICAL_MIRROR_PADDED_HILBERT"};

    for c = 1:2
        keys = double(trialKeys{c}(:));
        assert(numel(keys)==size(Env{c},2),'Trial-key count mismatch.');
        for t = 1:size(Env{c},2)
            amplitude = reshape(double(Env{c}(r,t,:)),1,[]);
            z = (amplitude-center)./scale;
            parents = findSegments(z>=thresholds(1));
            parentID = zeros(1,numel(z)); parentLeft = false(size(parents,1),1);
            for q = 1:size(parents,1)
                parentID(parents(q,1):parents(q,2)) = q;
                parentLeft(q) = parents(q,1)==1;
            end
            for h = 1:numel(thresholds)
                threshold = thresholds(h);
                segments = findSegments(z>=threshold);
                keep = (segments(:,2)-segments(:,1)+1)>=minSamples;
                segments = segments(keep,:);
                nAll = size(segments,1); nPrimary = 0;
                anyLeft = false; anyRight = false;
                for e = 1:nAll
                    first = segments(e,1); last = segments(e,2);
                    pid = parentID(first);
                    assert(pid>0,'Threshold event lacks z=1.5 parent.');
                    left = parentLeft(pid); right = last==numel(z); primary = ~left;
                    nPrimary = nPrimary+double(primary);
                    anyLeft = anyLeft||left; anyRight = anyRight||right;
                    zz = z(first:last); [peakZ,relativePeak] = max(zz);
                    peakSample = first+relativePeak-1; eventID = eventID+1;
                    eventRows(end+1,:) = {eventID,"CORTICAL",labels(r), ...
                        string(bandLabel),conditions(c),keys(t),t,threshold,pid, ...
                        first,last,peakSample,last-first+1,(last-first+1)/fs, ...
                        peakZ,sum(max(zz-threshold,0))/fs,left,right,primary}; %#ok<AGROW>
                end
                inventoryRows(end+1,:) = {"CORTICAL",labels(r),string(bandLabel), ...
                    conditions(c),keys(t),t,threshold,true,nAll>0,nPrimary>0, ...
                    nAll,nPrimary,anyLeft,anyRight,minSamples,minSamples/fs}; %#ok<AGROW>
            end
        end
    end
end
Events = cell2table(eventRows,'VariableNames',{'EventID','SignalType','Signal', ...
    'Band','Condition','OriginalTrial','TrialIndex','ZThreshold','ParentID', ...
    'OnsetSample','OffsetSample','PeakSample','DurationSamples','DurationSec', ...
    'PeakZ','AUCAboveThreshold','LeftCensored','RightCensored','PrimaryEligible'});
Inventory = cell2table(inventoryRows,'VariableNames',{'SignalType','Signal','Band', ...
    'Condition','OriginalTrial','TrialIndex','ZThreshold','TrialValid', ...
    'AnyBurstInclusive','AnyBurstPrimary','NBurstsInclusive','NBurstsPrimary', ...
    'AnyLeftCensored','AnyRightCensored','MinimumDurationSamples','MinimumDurationSec'});
Normalization = cell2table(normRows,'VariableNames',{'Signal','Band','NormCenter', ...
    'NormScale','NValues','NFinite','ScaleFallback','SignalType'});
end

function keys = getTrialKeys(P,conditions)
keys = cell(2,1); C = string(P.trialInventory.Condition);
for c = 1:2
    q = find(C==conditions(c));
    assert(numel(q)==size(P.roiData{c},2),'Trial inventory count mismatch.');
    keys{c} = double(P.trialInventory.OriginalTrial(q));
end
end

function Q = checkNesting(I)
key = {'Subject','Signal','Band','Condition','OriginalTrial'};
A = sortrows(I(I.ZThreshold==1.5,:),key);
B = sortrows(I(I.ZThreshold==2,:),key);
C = sortrows(I(I.ZThreshold==2.5,:),key);
assert(isequal(A(:,key),B(:,key)) && isequal(B(:,key),C(:,key)), ...
    'Nesting keys differ.');
Comparison = ["z2.0_without_z1.5";"z2.5_without_z2.0"];
NViolationsInclusive = [sum(B.AnyBurstInclusive&~A.AnyBurstInclusive); ...
    sum(C.AnyBurstInclusive&~B.AnyBurstInclusive)];
NViolationsPrimary = [sum(B.AnyBurstPrimary&~A.AnyBurstPrimary); ...
    sum(C.AnyBurstPrimary&~B.AnyBurstPrimary)];
NGroups = repmat(height(A),2,1);
Q = table(Comparison,NViolationsInclusive,NViolationsPrimary,NGroups);
end

function C = compareLegacyCounts(oldE,oldI,newE,newI)
bands = ["theta";"alpha";"beta"]; conditions = ["EMG2";"EMG3"];
rows = cell(18,10); row = 0;
for b = 1:3
    for c = 1:2
        for h = [1.5 2 2.5]
            row = row+1;
            oq = string(oldE.Band)==bands(b) & string(oldE.Condition)==conditions(c) & ...
                double(oldE.ZThreshold)==h & logical(oldE.PrimaryEligible);
            nq = string(newE.Band)==bands(b) & string(newE.Condition)==conditions(c) & ...
                double(newE.ZThreshold)==h & logical(newE.PrimaryEligible);
            oiq = string(oldI.Band)==bands(b) & string(oldI.Condition)==conditions(c) & ...
                double(oldI.ZThreshold)==h & logical(oldI.TrialValid);
            niq = string(newI.Band)==bands(b) & string(newI.Condition)==conditions(c) & ...
                double(newI.ZThreshold)==h & logical(newI.TrialValid);
            oldCount = sum(oq); newCount = sum(nq);
            rows(row,:) = {bands(b),conditions(c),h,sum(oiq),sum(niq), ...
                oldCount,newCount,newCount-oldCount,100*(newCount-oldCount)/max(oldCount,1), ...
                sum(niq)==sum(oiq)};
        end
    end
end
C = cell2table(rows,'VariableNames',{'Band','Condition','ZThreshold', ...
    'LegacyValidROITrials','MirrorValidROITrials','LegacyPrimaryEvents', ...
    'MirrorPrimaryEvents','EventCountDifference','RelativeDifferencePercent', ...
    'ValidTrialCountIdentical'});
assert(all(C.ValidTrialCountIdentical),'Valid trial counts changed unexpectedly.');
end

function B = loadBurstObject(filePath)
S = load(filePath); names = fieldnames(S); hits = false(numel(names),1);
for k = 1:numel(names)
    x = S.(names{k});
    hits(k) = isstruct(x) && isfield(x,'CorticalEvents') && ...
        isfield(x,'CorticalInventory');
end
assert(sum(hits)==1,'Expected exactly one burst object.');
B = S.(names{find(hits,1)});
end

function segments = findSegments(mask)
mask = logical(mask(:)'); d = diff([false mask false]);
segments = [find(d==1)' find(d==-1)'-1];
end

function envelope = analyticEnvelopeFFT(x)
% Absolute analytic signal using only base-MATLAB FFT/IFFT operations.
% This implements the same standard DFT-domain construction as hilbert(x)
% for a real vector while avoiding a Signal Processing Toolbox dependency.
x = double(x(:)');
n = numel(x);
assert(n>=2 && all(isfinite(x)),'Analytic-envelope input must be finite.');
multiplier = zeros(1,n);
if rem(n,2)==0
    multiplier(1) = 1;
    multiplier(n/2+1) = 1;
    multiplier(2:n/2) = 2;
else
    multiplier(1) = 1;
    multiplier(2:(n+1)/2) = 2;
end
analyticSignal = ifft(fft(x).*multiplier);
envelope = abs(analyticSignal);
end
