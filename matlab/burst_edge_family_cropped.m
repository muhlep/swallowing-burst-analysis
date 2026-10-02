function Unified = burst_edge_family_cropped( ...
    burstFile,legacyAllEdgesFile,varargin)
% Reconstruct participant-level co-burst effects and apply one max-|T|
% family across theta, alpha and beta (3 x 210 = 630 tests), separately
% for the mismatched-trial and non-zero circular-shift reference estimates.
%
% The function first requires exact reproduction of every locked group
% point estimate. If reproduction fails, it stops before interpreting the
% unified P values.

p = inputParser;
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'seed',20260918,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun,'Set confirmRun=true after checking both inputs.');

burstFile = char(string(burstFile));
legacyAllEdgesFile = char(string(legacyAllEdgesFile));
outputRoot = char(string(p.Results.outputRoot));
assert(isfile(burstFile),'Missing burst file: %s',burstFile);
assert(isfile(legacyAllEdgesFile),'Missing all-edge table: %s',legacyAllEdgesFile);
if isempty(outputRoot), outputRoot = fileparts(burstFile); end

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
bands = ["theta";"alpha";"beta"];
conditions = ["EMG2";"EMG3"];
assert(numel(subjects)==73,'Expected 73 healthy participants.');
assert(numel(rois)==21,'Expected 21 cortical ROIs.');
assert(all(ismember(bands,unique(string(CI.Band)))),'Standard bands missing.');
assert(all(ismember(conditions,unique(string(CI.Condition)))), ...
    'EMG2/EMG3 conditions missing.');

pairs = nchoosek(1:numel(rois),2);
nEdges = size(pairs,1);
nTests = numel(bands)*nEdges;
assert(nEdges==210 && nTests==630,'Unexpected family size.');
fs = 600;
nSamples = 600;

% Dtrial and Dcircular contain participant-level EMG2-minus-EMG3
% excess-occupancy differences in percentage points.
Dtrial = nan(numel(subjects),nTests);
Dcircular = nan(numel(subjects),nTests);
Band = strings(nTests,1); ROI_A = strings(nTests,1); ROI_B = strings(nTests,1);
test = 0;

for b = 1:numel(bands)
    band = bands(b);
    for e = 1:nEdges
        test = test+1;
        roiA = rois(pairs(e,1)); roiB = rois(pairs(e,2));
        Band(test) = band; ROI_A(test) = roiA; ROI_B(test) = roiB;
        for s = 1:numel(subjects)
            sid = subjects(s);
            excessTrial = nan(2,1);
            excessCircular = nan(2,1);
            for c = 1:2
                condition = conditions(c);
                [trialIDsA,masksA] = makeMasks(CI,CE,sid,roiA,band,condition,nSamples);
                [trialIDsB,masksB] = makeMasks(CI,CE,sid,roiB,band,condition,nSamples);
                common = intersect(trialIDsA,trialIDsB,'stable');
                assert(~isempty(common),'No common valid trials for %s %s %s.', ...
                    sid,band,condition);
                observed = zeros(numel(common),1);
                circular = zeros(numel(common),1);
                for k = 1:numel(common)
                    ia = find(trialIDsA==common(k),1);
                    ib = find(trialIDsB==common(k),1);
                    x = masksA(ia,:); y = masksB(ib,:);
                    overlap0 = sum(x & y)/nSamples;
                    occA = sum(x)/nSamples; occB = sum(y)/nSamples;
                    observed(k) = overlap0;
                    circular(k) = (nSamples*occA*occB-overlap0)/(nSamples-1);
                end
                obsMean = mean(observed);

                mismatchValues = zeros(0,1);
                for ia = 1:numel(trialIDsA)
                    keep = trialIDsB~=trialIDsA(ia);
                    if any(keep)
                        vals = mean(masksB(keep,:) & masksA(ia,:),2);
                        mismatchValues = [mismatchValues; vals]; %#ok<AGROW>
                    end
                end
                assert(~isempty(mismatchValues),'Mismatched-trial reference is empty.');
                excessTrial(c) = 100*(obsMean-mean(mismatchValues));
                excessCircular(c) = 100*(obsMean-mean(circular));
            end
            Dtrial(s,test) = excessTrial(1)-excessTrial(2);
            Dcircular(s,test) = excessCircular(1)-excessCircular(2);
        end
    end
end
assert(all(isfinite(Dtrial),'all') && all(isfinite(Dcircular),'all'), ...
    'Non-finite participant effects detected.');

% Exact point-estimate reproduction against the locked 630-edge table.
legacy = readtable(legacyAllEdgesFile,'TextType','string');
requiredLegacy = ["Band" "ROI_A" "ROI_B" ...
    "TrialSpecificDifferencePercent" "TimeLockedDifferencePercent"];
assert(all(ismember(requiredLegacy,string(legacy.Properties.VariableNames))), ...
    'Legacy table lacks required columns.');
assert(height(legacy)==630,'Expected 630 locked rows.');
legacyKey = edgeKey(legacy.Band,legacy.ROI_A,legacy.ROI_B);
newKey = edgeKey(Band,ROI_A,ROI_B);
[tf,loc] = ismember(newKey,legacyKey);
assert(all(tf) && numel(unique(loc))==630,'Locked edge keys do not match.');
trialDelta = mean(Dtrial,1)'-double(legacy.TrialSpecificDifferencePercent(loc));
circularDelta = mean(Dcircular,1)'-double(legacy.TimeLockedDifferencePercent(loc));
tolerance = 1e-10;
reproductionPass = all(abs(trialDelta)<tolerance) && ...
    all(abs(circularDelta)<tolerance);
assert(reproductionPass, ...
    ['Participant reconstruction does not reproduce the locked point estimates. ' ...
     'Do not interpret unified P values; inspect aggregation and validity masks.']);

rng(double(p.Results.seed),'twister');
[trialT,trialDz,trialRawP,trialMaxP] = maxT(Dtrial,double(p.Results.nPermutations));
rng(double(p.Results.seed),'twister');
[circT,circDz,circRawP,circMaxP] = maxT(Dcircular,double(p.Results.nPermutations));

trialDifference = mean(Dtrial,1)';
circularDifference = mean(Dcircular,1)';
trialPass = trialDifference>0 & isfinite(trialT) & trialMaxP<0.05;
circularPass = circularDifference>0 & isfinite(circT) & circMaxP<0.05;

UnifiedStats = table(Band,ROI_A,ROI_B,repmat(73,nTests,1), ...
    trialDifference,trialT,trialDz,trialRawP,trialMaxP,trialPass, ...
    circularDifference,circT,circDz,circRawP,circMaxP,circularPass, ...
    trialPass & circularPass, ...
    'VariableNames',{'Band','ROI_A','ROI_B','N', ...
    'TrialDifferencePercent','TrialT','TrialCohenDz','TrialRawP', ...
    'TrialUnifiedMaxTFWERP','TrialUnifiedMaxT05', ...
    'CircularDifferencePercent','CircularT','CircularCohenDz','CircularRawP', ...
    'CircularUnifiedMaxTFWERP','CircularUnifiedMaxT05', ...
    'BothUnifiedMaxT05'});

SubjectEffects = table;
for j = 1:nTests
    block = table(subjects,repmat(Band(j),73,1),repmat(ROI_A(j),73,1), ...
        repmat(ROI_B(j),73,1),Dtrial(:,j),Dcircular(:,j), ...
        'VariableNames',{'Subject','Band','ROI_A','ROI_B', ...
        'TrialDifferencePercent','CircularDifferencePercent'});
    SubjectEffects = [SubjectEffects; block]; %#ok<AGROW>
end

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
runID = ['BURST_EDGE_FAMILY_CROPPED_' stamp];
runDir = fullfile(outputRoot,runID); mkdir(runDir);
InputQC = table(73,21,3,210,630,2,fs,nSamples, ...
    double(p.Results.nPermutations),double(p.Results.seed),reproductionPass, ...
    max(abs(trialDelta)),max(abs(circularDelta)), ...
    'VariableNames',{'NSubjects','NROIs','NBands','NEdgesPerBand', ...
    'NTestsPerNull','ZThreshold','SamplingRateHz','WindowSamples', ...
    'NPermutations','Seed','LockedPointEstimatesReproduced', ...
    'MaximumTrialPointEstimateError','MaximumCircularPointEstimateError'});

writetable(InputQC,fullfile(runDir,'INPUT_QC.csv'));
writetable(UnifiedStats,fullfile(runDir,'UNIFIED_630_EDGE_STATS.csv'));
writetable(SubjectEffects,fullfile(runDir,'SUBJECT_LEVEL_630_EDGE_EFFECTS.csv'));
writetable(UnifiedStats(UnifiedStats.BothUnifiedMaxT05,:), ...
    fullfile(runDir,'EDGES_SURVIVING_BOTH_UNIFIED_CORRECTIONS.csv'));

Unified = struct;
Unified.cfg = struct('runID',runID,'runDir',runDir,'burstFile',burstFile, ...
    'legacyAllEdgesFile',legacyAllEdgesFile,'created',stamp, ...
    'nPermutations',double(p.Results.nPermutations),'seed',double(p.Results.seed));
Unified.InputQC = InputQC;
Unified.UnifiedStats = UnifiedStats;
Unified.SubjectEffects = SubjectEffects;
save(fullfile(runDir,[runID '_RESULTS.mat']),'Unified','-v7.3');
end

function [trialIDs,masks] = makeMasks(CI,CE,sid,roi,band,condition,nSamples)
q = string(CI.Subject)==sid & string(CI.Signal)==roi & ...
    string(CI.Band)==band & string(CI.Condition)==condition & ...
    double(CI.ZThreshold)==2;
if ismember('TrialValid',CI.Properties.VariableNames)
    q = q & logical(CI.TrialValid);
end
trialIDs = unique(double(CI.OriginalTrial(q)),'stable');
assert(~isempty(trialIDs),'No valid trials for %s %s %s %s.',sid,roi,band,condition);
masks = false(numel(trialIDs),nSamples);
e = string(CE.Subject)==sid & string(CE.Signal)==roi & ...
    string(CE.Band)==band & string(CE.Condition)==condition & ...
    double(CE.ZThreshold)==2 & logical(CE.PrimaryEligible);
E = CE(e,:);
for k = 1:height(E)
    row = find(trialIDs==double(E.OriginalTrial(k)),1);
    if isempty(row), continue, end
    first = max(1,round(double(E.OnsetSample(k))));
    last = min(nSamples,round(double(E.OffsetSample(k))));
    if first<=last, masks(row,first:last) = true; end
end
end

function B = loadBurstObject(filePath)
S = load(filePath);
names = fieldnames(S);
hits = false(numel(names),1);
for k = 1:numel(names)
    x = S.(names{k});
    hits(k) = isstruct(x) && isfield(x,'CorticalEvents') && ...
        isfield(x,'CorticalInventory');
end
assert(sum(hits)==1,'Expected exactly one burst object with events/inventory.');
B = S.(names{find(hits,1)});
end

function key = edgeKey(band,a,b)
band = string(band); a = string(a); b = string(b);
lo = strings(size(a)); hi = strings(size(a));
for k = 1:numel(a)
    pair = sort([a(k) b(k)]);
    lo(k) = pair(1); hi(k) = pair(2);
end
key = band + "|" + lo + "|" + hi;
end

function [obsT,dz,rawP,maxP] = maxT(D,nPerm)
[n,m] = size(D);
mu = mean(D,1);
sd = std(D,0,1);
valid = isfinite(mu) & isfinite(sd) & sd>0;
obsT = zeros(m,1);
dz = zeros(m,1);
obsT(valid) = (mu(valid)./(sd(valid)/sqrt(n)))';
dz(valid) = (mu(valid)./sd(valid))';
rawCount = zeros(1,m);
maxCount = zeros(1,m);
sumSq = sum(D.^2,1);
batchSize = 1000;
done = 0;
while done<nPerm
    b = min(batchSize,nPerm-done);
    signs = 2*(rand(b,n)>0.5)-1;
    permMean = (signs*D)/n;
    numerator = max(sumSq-n*(permMean.^2),0);
    permSD = sqrt(numerator/(n-1));
    permT = zeros(b,m);
    permT(:,valid) = permMean(:,valid)./(permSD(:,valid)/sqrt(n));
    rawCount = rawCount + sum(abs(permT)>=abs(obsT'),1);
    mx = max(abs(permT),[],2);
    maxCount = maxCount + sum(mx>=abs(obsT'),1);
    done = done+b;
end
rawP = ((rawCount+1)/(nPerm+1))';
maxP = ((maxCount+1)/(nPerm+1))';
rawP(~valid) = 1;
maxP(~valid) = 1;
end
