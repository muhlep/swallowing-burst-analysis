function Audit = burst_bounded_onset_lag( ...
    burstFile,finalEdgeSource,varargin)
% Bounded onset-lag sensitivity for the final six beta co-burst edges.
% It does not estimate direction, propagation, hierarchy or causality.
%
% Required inputs
%   burstFile       MAT file containing a struct with CorticalEvents and
%                   CorticalInventory (for example BurstFullR02).
%   finalEdgeSource MAT or CSV containing the final corrected edge table.
%                   Only beta edges passing both null models are retained.
%
% Name-value inputs
%   outputRoot      Parent directory for a new versioned result folder.
%   legacyStatsFile Optional legacy ONSETLAGSTATS CSV for numerical audit.
%   nPermutations   Sign-flip randomizations (default 100000).
%   nBootstrap      Participant bootstrap samples (default 20000).
%   seed            Reproducible random seed (default 20260913).
%   confirmRun      Must be true.

p = inputParser;
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'legacyStatsFile','',@(x)ischar(x)||isstring(x));
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
addParameter(p,'seed',20260913,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun, ...
    'Set confirmRun=true after checking burstFile and finalEdgeSource.');

burstFile = char(string(burstFile));
finalEdgeSource = char(string(finalEdgeSource));
legacyStatsFile = char(string(p.Results.legacyStatsFile));
assert(isfile(burstFile),'Burst input is missing: %s',burstFile);
assert(isfile(finalEdgeSource),'Final edge source is missing: %s',finalEdgeSource);
if ~isempty(legacyStatsFile)
    assert(isfile(legacyStatsFile),'Legacy stats file is missing: %s',legacyStatsFile);
end

%% Load and lock the healthy event-level input
B = loadBurstObject(burstFile);
CE = B.CorticalEvents;
CI = B.CorticalInventory;
requiredEvent = ["Subject" "Signal" "Band" "Condition" "ZThreshold" ...
    "OriginalTrial" "OnsetSample" "OffsetSample" "PrimaryEligible"];
requiredInventory = ["Subject" "Signal" "Band" "Condition" ...
    "ZThreshold" "OriginalTrial"];
assert(all(ismember(requiredEvent,string(CE.Properties.VariableNames))), ...
    'CorticalEvents lacks required columns.');
assert(all(ismember(requiredInventory,string(CI.Properties.VariableNames))), ...
    'CorticalInventory lacks required columns.');

subjects = sort(unique(string(CI.Subject)));
rois = sort(unique(string(CI.Signal)));
assert(numel(subjects)==73,'Expected 73 healthy participants, found %d.',numel(subjects));
assert(numel(rois)==21,'Expected 21 cortical ROIs, found %d.',numel(rois));
assert(all(~ismissing(subjects) & strlength(subjects)>0), ...
    'Participant pairing keys must be non-missing and non-empty.');
conditions = ["EMG2";"EMG3"];
assert(all(ismember(conditions,unique(string(CI.Condition)))), ...
    'EMG2/EMG3 condition pair is incomplete.');

%% Load the final six beta edges and verify their identity
edges = loadFinalBetaEdges(finalEdgeSource);
expected = [ ...
    "BA44_R"    "Ins_Ant_R"; ...
    "Ins_Ant_L" "Ins_Post_L"; ...
    "Ins_Ant_R" "Ins_Post_R"; ...
    "Ins_Post_R" "PMd_R"; ...
    "Ins_Post_R" "S2_R"; ...
    "Ins_Post_R" "SMG_R"];
assert(height(edges)==6,'Expected exactly six final beta edges, found %d.',height(edges));
assert(isequal(sort(edgeKeys(edges.ROI_A,edges.ROI_B)), ...
    sort(edgeKeys(expected(:,1),expected(:,2)))), ...
    'Final beta-edge identity differs from the locked six-edge signature.');
assert(all(ismember(string(edges.ROI_A),rois)) && ...
    all(ismember(string(edges.ROI_B),rois)), ...
    'At least one final edge contains an unknown ROI.');

%% Reconstruct one-to-one overlapping-event matches
fs = 600;
nearSamples = round(0.010*fs);
delayedMaxSamples = round(0.100*fs);
nSubjects = numel(subjects);
nEdges = height(edges);
nRows = nSubjects*numel(conditions)*nEdges;
rows = cell(nRows,15);
row = 0;

for s = 1:nSubjects
    sid = subjects(s);
    for c = 1:numel(conditions)
        condition = conditions(c);
        iq = string(CI.Subject)==sid & string(CI.Band)=="beta" & ...
            string(CI.Condition)==condition & double(CI.ZThreshold)==2;
        trials = unique(double(CI.OriginalTrial(iq)));
        assert(~isempty(trials),'No beta z=2 trials for %s %s.',sid,condition);
        eq = string(CE.Subject)==sid & string(CE.Band)=="beta" & ...
            string(CE.Condition)==condition & double(CE.ZThreshold)==2 & ...
            logical(CE.PrimaryEligible);
        events = CE(eq,:);

        for e = 1:nEdges
            roiA = string(edges.ROI_A(e));
            roiB = string(edges.ROI_B(e));
            nMatched = 0;
            nNear = 0;
            nDelayed = 0;
            nLong = 0;
            absoluteLags = zeros(0,1);
            overlaps = zeros(0,1);

            for t = 1:numel(trials)
                trial = trials(t);
                A = events(string(events.Signal)==roiA & ...
                    double(events.OriginalTrial)==trial,:);
                Btrial = events(string(events.Signal)==roiB & ...
                    double(events.OriginalTrial)==trial,:);
                pairs = matchOverlappingBursts(A,Btrial);
                for k = 1:size(pairs,1)
                    lag = abs(double(Btrial.OnsetSample(pairs(k,2))) - ...
                        double(A.OnsetSample(pairs(k,1))));
                    nMatched = nMatched+1;
                    nNear = nNear+(lag<=nearSamples);
                    nDelayed = nDelayed+(lag>nearSamples && lag<=delayedMaxSamples);
                    nLong = nLong+(lag>delayedMaxSamples);
                    absoluteLags(end+1,1) = lag; %#ok<AGROW>
                    overlaps(end+1,1) = pairs(k,3); %#ok<AGROW>
                end
            end

            if isempty(absoluteLags)
                medianAbsoluteLagSec = NaN;
                medianOverlapSec = NaN;
            else
                medianAbsoluteLagSec = median(absoluteLags)/fs;
                medianOverlapSec = median(overlaps)/fs;
            end
            row = row+1;
            rows(row,:) = {sid,condition,roiA,roiB,numel(trials), ...
                nMatched,nNear,nDelayed,nLong,nNear/numel(trials), ...
                nDelayed/numel(trials),nLong/numel(trials), ...
                medianAbsoluteLagSec,medianOverlapSec,2};
        end
    end
end

SubjectEdgeLag = cell2table(rows,'VariableNames',{ ...
    'Subject','Condition','ROI_A','ROI_B','NTrials','NMatched', ...
    'NNearZero','NDelayed','NLong','NearZeroMatchesPerTrial', ...
    'DelayedMatchesPerTrial','LongMatchesPerTrial', ...
    'MedianAbsoluteLagSec','MedianOverlapSec','ZThreshold'});

%% Participant-level EMG2-minus-EMG3 contrasts
Dnear = nan(nSubjects,nEdges);
Ddelayed = nan(nSubjects,nEdges);
Dlong = nan(nSubjects,nEdges);
nMatchedEMG2 = nan(nEdges,1);
medianAbsoluteLagEMG2 = nan(nEdges,1);

for e = 1:nEdges
    q = edgeMatch(SubjectEdgeLag.ROI_A,SubjectEdgeLag.ROI_B, ...
        edges.ROI_A(e),edges.ROI_B(e));
    A = sortrows(SubjectEdgeLag(q & ...
        string(SubjectEdgeLag.Condition)=="EMG2",:),'Subject');
    R = sortrows(SubjectEdgeLag(q & ...
        string(SubjectEdgeLag.Condition)=="EMG3",:),'Subject');
    assert(height(A)==73 && height(R)==73 && ...
        isequal(string(A.Subject),subjects) && ...
        isequal(string(R.Subject),subjects), ...
        'Participant pairing failed for edge %s--%s.', ...
        string(edges.ROI_A(e)),string(edges.ROI_B(e)));
    Dnear(:,e) = 100*(A.NearZeroMatchesPerTrial-R.NearZeroMatchesPerTrial);
    Ddelayed(:,e) = 100*(A.DelayedMatchesPerTrial-R.DelayedMatchesPerTrial);
    Dlong(:,e) = 100*(A.LongMatchesPerTrial-R.LongMatchesPerTrial);
    nMatchedEMG2(e) = sum(A.NMatched);
    pooled = SubjectEdgeLag(q & ...
        string(SubjectEdgeLag.Condition)=="EMG2",:);
    weights = double(pooled.NMatched);
    valid = weights>0 & isfinite(pooled.MedianAbsoluteLagSec);
    if any(valid)
        medianAbsoluteLagEMG2(e) = weightedMedian( ...
            double(pooled.MedianAbsoluteLagSec(valid)),weights(valid));
    end
end

% The delayed 10-100 ms component is the prespecified primary leakage
% sensitivity family. Near-zero and >100 ms components are secondary
% diagnostics, each corrected across the same six final edges.
rng(double(p.Results.seed),'twister');
[delayedStats,delayedCI] = maxTInference(Ddelayed, ...
    double(p.Results.nPermutations),double(p.Results.nBootstrap));
[nearStats,nearCI] = maxTInference(Dnear, ...
    double(p.Results.nPermutations),double(p.Results.nBootstrap));
[longStats,longCI] = maxTInference(Dlong, ...
    double(p.Results.nPermutations),double(p.Results.nBootstrap));

EdgeLagStats = table(string(edges.ROI_A),string(edges.ROI_B), ...
    repmat(73,nEdges,1),nMatchedEMG2,medianAbsoluteLagEMG2, ...
    mean(Ddelayed,1)',delayedCI(:,1),delayedCI(:,2),delayedStats.Dz, ...
    delayedStats.RawP,delayedStats.MaxTP,delayedStats.MaxTP<0.05, ...
    mean(Dnear,1)',nearCI(:,1),nearCI(:,2),nearStats.Dz, ...
    nearStats.RawP,nearStats.MaxTP,nearStats.MaxTP<0.05, ...
    mean(Dlong,1)',longCI(:,1),longCI(:,2),longStats.Dz, ...
    longStats.RawP,longStats.MaxTP,longStats.MaxTP<0.05, ...
    'VariableNames',{ ...
    'ROI_A','ROI_B','N','NMatchedEventsEMG2', ...
    'WeightedMedianParticipantMedianAbsoluteLagEMG2Sec', ...
    'Delayed10to100msDifferencePercent','DelayedCI95Lower', ...
    'DelayedCI95Upper','DelayedCohenDz','DelayedRawP', ...
    'DelayedMaxTFWERP','DelayedMaxT05', ...
    'NearZeroDifferencePercent','NearZeroCI95Lower','NearZeroCI95Upper', ...
    'NearZeroCohenDz','NearZeroRawP','NearZeroMaxTFWERP','NearZeroMaxT05', ...
    'LongerThan100msDifferencePercent','LongCI95Lower','LongCI95Upper', ...
    'LongCohenDz','LongRawP','LongMaxTFWERP','LongMaxT05'});

EdgeLagStats = sortrows(EdgeLagStats,{'ROI_A','ROI_B'});

%% Optional exact comparison with the legacy seven-edge output
LegacyComparison = table();
legacyPass = true;
if ~isempty(legacyStatsFile)
    old = readtable(legacyStatsFile,'TextType','string');
    requiredLegacy = ["ROI_A" "ROI_B" "NMatchedEventsEMG2" ...
        "EMG2MinusEMG3NearZeroRatePercent" ...
        "EMG2MinusEMG3PhysiologicLagRatePercent"];
    assert(all(ismember(requiredLegacy,string(old.Properties.VariableNames))), ...
        'Legacy stats file lacks required columns.');
    compareRows = cell(nEdges,8);
    for e = 1:nEdges
        qOld = edgeMatch(old.ROI_A,old.ROI_B, ...
            EdgeLagStats.ROI_A(e),EdgeLagStats.ROI_B(e));
        assert(sum(qOld)==1,'Legacy row missing or duplicated for %s--%s.', ...
            EdgeLagStats.ROI_A(e),EdgeLagStats.ROI_B(e));
        o = old(qOld,:);
        countMatch = double(o.NMatchedEventsEMG2)== ...
            EdgeLagStats.NMatchedEventsEMG2(e);
        nearDelta = EdgeLagStats.NearZeroDifferencePercent(e)- ...
            double(o.EMG2MinusEMG3NearZeroRatePercent);
        delayedDelta = EdgeLagStats.Delayed10to100msDifferencePercent(e)- ...
            double(o.EMG2MinusEMG3PhysiologicLagRatePercent);
        effectMatch = abs(nearDelta)<1e-10 && abs(delayedDelta)<1e-10;
        compareRows(e,:) = {EdgeLagStats.ROI_A(e),EdgeLagStats.ROI_B(e), ...
            countMatch,nearDelta,delayedDelta,effectMatch, ...
            EdgeLagStats.DelayedMaxT05(e),EdgeLagStats.NearZeroMaxT05(e)};
    end
    LegacyComparison = cell2table(compareRows,'VariableNames',{ ...
        'ROI_A','ROI_B','MatchedEventCountExact','NearEffectDifference', ...
        'DelayedEffectDifference','LegacyEffectExact', ...
        'DelayedMaxT05SixEdgeFamily','NearZeroMaxT05SixEdgeFamily'});
    legacyPass = all(LegacyComparison.MatchedEventCountExact) && ...
        all(LegacyComparison.LegacyEffectExact);
end

%% Audit decision and versioned output
InputQC = table(73,21,6,2,fs,nearSamples,delayedMaxSamples, ...
    numel(unique(subjects))==73,legacyPass, ...
    'VariableNames',{'NSubjects','NROIs','NFinalBetaEdges','ZThreshold', ...
    'SamplingRateHz','NearZeroMaximumSamples','DelayedMaximumSamples', ...
    'UniqueOpaqueParticipantKeys','LegacyReproductionPass'});

ClaimBoundary = table( ...
    ["ALLOWED";"ALLOWED_IF_REPRODUCED";"FORBIDDEN";"FORBIDDEN"], ...
    ["Lag-class sensitivity was restricted to the final six beta edges."; ...
     "A swallowing-related delayed-overlap component was present for the specified edges."; ...
     "The signed lag identifies directed information flow or anatomical propagation."; ...
     "The result establishes a cortical hierarchy or pacemaker."], ...
    'VariableNames',{'Status','Wording'});

outputRoot = char(string(p.Results.outputRoot));
if isempty(outputRoot)
    outputRoot = fileparts(burstFile);
end
timeTag = char(datetime('now','Format','yyyyMMdd_HHmmss'));
runID = ['BURST_BOUNDED_ONSET_LAG_' timeTag];
runDir = fullfile(outputRoot,runID);
assert(~isfolder(runDir),'Refusing to overwrite existing output: %s',runDir);
[ok,msg] = mkdir(runDir);
assert(ok,'Cannot create output directory: %s',msg);

cfg = struct('runID',runID,'runDir',string(runDir), ...
    'burstFile',string(burstFile),'finalEdgeSource',string(finalEdgeSource), ...
    'legacyStatsFile',string(legacyStatsFile),'band','beta', ...
    'zThreshold',2,'fs',fs,'nearZeroDefinition','absolute onset lag <=10 ms', ...
    'delayedDefinition','10 < absolute onset lag <=100 ms', ...
    'matching','greedy one-to-one maximum-overlap matching within trial', ...
    'primaryFamily','delayed component; max-|T| FWER across six final beta edges', ...
    'secondaryFamilies',['near-zero and >100-ms diagnostics; separate ' ...
        'max-|T| FWER across the same six edges'], ...
    'nPermutations',p.Results.nPermutations, ...
    'nBootstrap',p.Results.nBootstrap,'seed',p.Results.seed, ...
    'interpretation',['spatial-leakage sensitivity only; no direction, ' ...
        'propagation, hierarchy or causality']);

Audit = struct('InputQC',InputQC,'FinalEdges',edges, ...
    'SubjectEdgeLag',SubjectEdgeLag,'EdgeLagStats',EdgeLagStats, ...
    'LegacyComparison',LegacyComparison,'ClaimBoundary',ClaimBoundary, ...
    'cfg',cfg);
save(fullfile(runDir,[runID '_AUDIT.mat']),'Audit','-v7.3');
writetable(InputQC,fullfile(runDir,'INPUT_QC.csv'));
writetable(edges,fullfile(runDir,'FINAL_SIX_BETA_EDGES.csv'));
writetable(SubjectEdgeLag,fullfile(runDir,'SUBJECT_EDGE_LAG.csv'));
writetable(EdgeLagStats,fullfile(runDir,'EDGE_LAG_STATS.csv'));
writetable(ClaimBoundary,fullfile(runDir,'CLAIM_BOUNDARY.csv'));
if ~isempty(LegacyComparison)
    writetable(LegacyComparison,fullfile(runDir,'LEGACY_REPRODUCTION.csv'));
end

fprintf('\nBounded onset-lag analysis complete:\n%s\n',runDir);
disp(InputQC)
fprintf('\nFINAL SIX-EDGE LAG-CLASS RESULTS\n')
disp(EdgeLagStats)
fprintf('\nDelayed component supported after six-edge max-T: %d/6\n', ...
    sum(EdgeLagStats.DelayedMaxT05));
fprintf(['Interpretation guard: this is a leakage sensitivity only. ' ...
    'Do not infer direction, propagation, hierarchy or causality.\n']);
end

function B = loadBurstObject(file)
S = load(file);
names = fieldnames(S);
hits = false(numel(names),1);
for k = 1:numel(names)
    x = S.(names{k});
    hits(k) = isstruct(x) && isfield(x,'CorticalEvents') && ...
        isfield(x,'CorticalInventory') && istable(x.CorticalEvents) && ...
        istable(x.CorticalInventory);
end
assert(sum(hits)==1, ...
    'Expected exactly one burst object with CorticalEvents/CorticalInventory.');
B = S.(names{find(hits,1)});
end

function edges = loadFinalBetaEdges(file)
[~,~,ext] = fileparts(file);
if strcmpi(ext,'.csv')
    T = readtable(file,'TextType','string');
else
    S = load(file);
    T = findEdgeTable(S);
end
vars = string(T.Properties.VariableNames);
assert(all(ismember(["ROI_A" "ROI_B"],vars)), ...
    'Final edge source lacks ROI_A/ROI_B.');
q = true(height(T),1);
if ismember("Band",vars)
    q = q & lower(string(T.Band))=="beta";
end
if ismember("BothNullsPass",vars)
    q = q & logical(T.BothNullsPass);
elseif ismember("BothNullDifferencesMaxT05",vars)
    q = q & logical(T.BothNullDifferencesMaxT05);
else
    error('Final edge table lacks a both-null pass flag.');
end
edges = T(q,{'ROI_A','ROI_B'});
edges.ROI_A = string(edges.ROI_A);
edges.ROI_B = string(edges.ROI_B);
edges = sortrows(edges,{'ROI_A','ROI_B'});
end

function T = findEdgeTable(S)
T = table();
names = fieldnames(S);
for k = 1:numel(names)
    x = S.(names{k});
    if istable(x) && all(ismember({'ROI_A','ROI_B'},x.Properties.VariableNames))
        T = x;
        return
    end
    if isstruct(x)
        sub = fieldnames(x);
        for j = 1:numel(sub)
            y = x.(sub{j});
            if istable(y) && all(ismember({'ROI_A','ROI_B'}, ...
                    y.Properties.VariableNames)) && ...
                    (ismember('BothNullsPass',y.Properties.VariableNames) || ...
                     ismember('BothNullDifferencesMaxT05',y.Properties.VariableNames))
                T = y;
                return
            end
        end
    end
end
error('No corrected edge table was found in the MAT file.');
end

function keys = edgeKeys(A,B)
A = string(A); B = string(B);
assert(isequal(size(A),size(B)),'Edge endpoint arrays differ in size.');
keys = strings(size(A));
for i = 1:numel(A)
    pair = sort([A(i) B(i)]);
    keys(i) = pair(1)+"|"+pair(2);
end
end

function q = edgeMatch(A,B,a,b)
q = edgeKeys(A,B)==edgeKeys(a,b);
end

function pairs = matchOverlappingBursts(A,B)
pairs = zeros(0,3);
if isempty(A) || isempty(B)
    return
end
candidates = zeros(0,4);
for i = 1:height(A)
    for j = 1:height(B)
        overlap = min(double(A.OffsetSample(i)),double(B.OffsetSample(j))) - ...
            max(double(A.OnsetSample(i)),double(B.OnsetSample(j)))+1;
        if overlap>0
            lag = abs(double(B.OnsetSample(j))-double(A.OnsetSample(i)));
            candidates(end+1,:) = [i j overlap lag]; %#ok<AGROW>
        end
    end
end
if isempty(candidates)
    return
end
candidates = sortrows(candidates,[-3 4]);
usedA = false(height(A),1);
usedB = false(height(B),1);
for k = 1:size(candidates,1)
    i = candidates(k,1);
    j = candidates(k,2);
    if ~usedA(i) && ~usedB(j)
        pairs(end+1,:) = candidates(k,1:3); %#ok<AGROW>
        usedA(i) = true;
        usedB(j) = true;
    end
end
end

function [O,CI] = maxTInference(D,nPerm,nBoot)
n = size(D,1);
k = size(D,2);
assert(all(isfinite(D),'all'),'Non-finite participant contrasts detected.');
mu = mean(D,1);
sd = std(D,0,1);
tObs = zeros(1,k);
valid = sd>0 & isfinite(sd);
tObs(valid) = mu(valid)./(sd(valid)/sqrt(n));
ss = sum(D.^2,1);
raw = zeros(1,k);
adjusted = zeros(1,k);
chunk = 1000;
for first = 1:chunk:nPerm
    m = min(chunk,nPerm-first+1);
    signs = 2*(rand(m,n)>0.5)-1;
    nullMean = (signs*D)/n;
    nullVar = (repmat(ss,m,1)-n*nullMean.^2)/(n-1);
    nullT = zeros(size(nullMean));
    ok = nullVar>0;
    nullT(ok) = nullMean(ok)./sqrt(nullVar(ok)/n);
    maxAbsT = max(abs(nullT),[],2);
    raw = raw+sum(abs(nullT)>=abs(tObs),1);
    adjusted = adjusted+sum(maxAbsT>=abs(tObs),1);
end
O = struct('Dz',(mu./sd)', ...
    'RawP',((raw+1)/(nPerm+1))', ...
    'MaxTP',((adjusted+1)/(nPerm+1))');
O.Dz(~isfinite(O.Dz)) = 0;

CI = nan(k,2);
for j = 1:k
    d = D(:,j);
    boot = zeros(nBoot,1);
    done = 0;
    while done<nBoot
        m = min(2000,nBoot-done);
        idx = randi(n,n,m);
        sampled = reshape(d(idx),n,m);
        boot(done+(1:m)) = mean(sampled,1)';
        done = done+m;
    end
    CI(j,:) = prctile(boot,[2.5 97.5]);
end
end

function value = weightedMedian(x,w)
[x,order] = sort(x);
w = w(order);
cut = 0.5*sum(w);
value = x(find(cumsum(w)>=cut,1,'first'));
end
