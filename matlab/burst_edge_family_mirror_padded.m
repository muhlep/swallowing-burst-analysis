function Results = burst_edge_family_mirror_padded( ...
    mirrorBurstFile,legacyUnifiedFile,varargin)
% Recompute the complete theta/alpha/beta co-burst edge family from the
% mirror-padded analytic-amplitude event masks and compare them with the
% cropped-envelope 630-edge family.

p = inputParser;
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'seed',20260924,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun, ...
    'Set confirmRun=true only after checking both input files.');

mirrorBurstFile = char(string(mirrorBurstFile));
legacyUnifiedFile = char(string(legacyUnifiedFile));
outputRoot = char(string(p.Results.outputRoot));
assert(isfile(mirrorBurstFile),'Missing mirror-padded burst file: %s', ...
    mirrorBurstFile);
assert(isfile(legacyUnifiedFile),'Missing cropped-envelope edge-family file: %s', ...
    legacyUnifiedFile);
if isempty(outputRoot), outputRoot = fileparts(mirrorBurstFile); end
assert(isfolder(outputRoot),'Output root does not exist: %s',outputRoot);

B = loadBurstObject(mirrorBurstFile);
CE = B.CorticalEvents;
CI = B.CorticalInventory;
Legacy = loadUnifiedObject(legacyUnifiedFile);

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
assert(all(ismember(bands,unique(string(CI.Band)))), ...
    'One or more standard bands are missing.');
assert(all(ismember(conditions,unique(string(CI.Condition)))), ...
    'EMG2/EMG3 conditions are missing.');

pairs = nchoosek(1:numel(rois),2);
nEdges = size(pairs,1);
nTests = numel(bands)*nEdges;
nSubjects = numel(subjects);
nSamples = inferWindowSamples(CI,CE);
fs = 600;
assert(nEdges==210 && nTests==630,'Unexpected family size.');
assert(nSamples==600,'Expected a 600-sample analysis window.');

% Participant-level values are EMG2-minus-EMG3 differences in excess
% co-burst occupancy, expressed in percentage points.
Dtrial = nan(nSubjects,nTests);
Dcircular = nan(nSubjects,nTests);
Band = strings(nTests,1);
ROI_A = strings(nTests,1);
ROI_B = strings(nTests,1);
test = 0;

for b = 1:numel(bands)
    band = bands(b);
    fprintf('\nMirror-padded edge family: %s\n',upper(band));
    for e = 1:nEdges
        test = test+1;
        roiA = rois(pairs(e,1));
        roiB = rois(pairs(e,2));
        Band(test) = band;
        ROI_A(test) = roiA;
        ROI_B(test) = roiB;
        if e==1 || mod(e,25)==0 || e==nEdges
            fprintf('  edge %d/%d: %s - %s\n',e,nEdges,roiA,roiB);
        end
        for s = 1:nSubjects
            sid = subjects(s);
            excessTrial = nan(2,1);
            excessCircular = nan(2,1);
            for c = 1:2
                condition = conditions(c);
                [trialIDsA,masksA] = makeMasks( ...
                    CI,CE,sid,roiA,band,condition,nSamples);
                [trialIDsB,masksB] = makeMasks( ...
                    CI,CE,sid,roiB,band,condition,nSamples);
                common = intersect(trialIDsA,trialIDsB,'stable');
                assert(~isempty(common), ...
                    'No common valid trials for %s %s %s.', ...
                    sid,band,condition);

                observed = zeros(numel(common),1);
                circular = zeros(numel(common),1);
                for k = 1:numel(common)
                    ia = find(trialIDsA==common(k),1);
                    ib = find(trialIDsB==common(k),1);
                    x = masksA(ia,:);
                    y = masksB(ib,:);
                    overlap0 = sum(x & y)/nSamples;
                    occA = sum(x)/nSamples;
                    occB = sum(y)/nSamples;
                    observed(k) = overlap0;
                    circular(k) = ...
                        (nSamples*occA*occB-overlap0)/(nSamples-1);
                end
                obsMean = mean(observed);

                mismatchSum = 0;
                mismatchN = 0;
                for ia = 1:numel(trialIDsA)
                    keep = trialIDsB~=trialIDsA(ia);
                    if any(keep)
                        values = mean(masksB(keep,:) & masksA(ia,:),2);
                        mismatchSum = mismatchSum + sum(values);
                        mismatchN = mismatchN + numel(values);
                    end
                end
                assert(mismatchN>0,'Mismatched-trial reference is empty.');
                mismatchMean = mismatchSum/mismatchN;
                excessTrial(c) = 100*(obsMean-mismatchMean);
                excessCircular(c) = 100*(obsMean-mean(circular));
            end
            Dtrial(s,test) = excessTrial(1)-excessTrial(2);
            Dcircular(s,test) = excessCircular(1)-excessCircular(2);
        end
    end
end
assert(all(isfinite(Dtrial),'all') && all(isfinite(Dcircular),'all'), ...
    'Non-finite participant effects detected.');

rng(double(p.Results.seed),'twister');
[trialT,trialDz,trialRawP,trialMaxP] = ...
    maxT(Dtrial,double(p.Results.nPermutations));
rng(double(p.Results.seed),'twister');
[circT,circDz,circRawP,circMaxP] = ...
    maxT(Dcircular,double(p.Results.nPermutations));

trialDifference = mean(Dtrial,1)';
circularDifference = mean(Dcircular,1)';
trialPass = trialDifference>0 & isfinite(trialT) & trialMaxP<0.05;
circularPass = circularDifference>0 & isfinite(circT) & circMaxP<0.05;

MirrorStats = table(Band,ROI_A,ROI_B,repmat(nSubjects,nTests,1), ...
    trialDifference,trialT,trialDz,trialRawP,trialMaxP,trialPass, ...
    circularDifference,circT,circDz,circRawP,circMaxP,circularPass, ...
    trialPass & circularPass, ...
    'VariableNames',{'Band','ROI_A','ROI_B','N', ...
    'TrialDifferencePercent','TrialT','TrialCohenDz','TrialRawP', ...
    'TrialUnifiedMaxTFWERP','TrialUnifiedMaxT05', ...
    'CircularDifferencePercent','CircularT','CircularCohenDz', ...
    'CircularRawP','CircularUnifiedMaxTFWERP','CircularUnifiedMaxT05', ...
    'BothUnifiedMaxT05'});

SubjectEffects = table;
for j = 1:nTests
    block = table(subjects,repmat(Band(j),nSubjects,1), ...
        repmat(ROI_A(j),nSubjects,1),repmat(ROI_B(j),nSubjects,1), ...
        Dtrial(:,j),Dcircular(:,j), ...
        'VariableNames',{'Subject','Band','ROI_A','ROI_B', ...
        'TrialDifferencePercent','CircularDifferencePercent'});
    SubjectEffects = [SubjectEffects; block]; %#ok<AGROW>
end

% Align with the cropped-envelope family and classify candidate-set changes.
LegacyStats = Legacy.UnifiedStats;
assert(height(LegacyStats)==630,'Cropped-envelope table must contain 630 rows.');
requiredLegacy = ["Band" "ROI_A" "ROI_B" ...
    "TrialDifferencePercent" "TrialUnifiedMaxTFWERP" ...
    "CircularDifferencePercent" "CircularUnifiedMaxTFWERP" ...
    "BothUnifiedMaxT05"];
assert(all(ismember(requiredLegacy, ...
    string(LegacyStats.Properties.VariableNames))), ...
    'Cropped-envelope table lacks required columns.');
legacyKey = edgeKey(LegacyStats.Band,LegacyStats.ROI_A,LegacyStats.ROI_B);
mirrorKey = edgeKey(MirrorStats.Band,MirrorStats.ROI_A,MirrorStats.ROI_B);
[tf,loc] = ismember(mirrorKey,legacyKey);
assert(all(tf) && numel(unique(loc))==630, ...
    'Locked and mirror-padded edge keys do not match.');
LegacyStats = LegacyStats(loc,:);

legacyBoth = logical(LegacyStats.BothUnifiedMaxT05);
mirrorBoth = logical(MirrorStats.BothUnifiedMaxT05);
status = repmat("NONCANDIDATE_BOTH",nTests,1);
status(legacyBoth & mirrorBoth) = "RETAINED_AFTER_MIRROR_PADDING";
status(legacyBoth & ~mirrorBoth) = "LOST_AFTER_MIRROR_PADDING";
status(~legacyBoth & mirrorBoth) = "GAINED_AFTER_MIRROR_PADDING";

Comparison = table(Band,ROI_A,ROI_B, ...
    double(LegacyStats.TrialDifferencePercent),trialDifference, ...
    trialDifference-double(LegacyStats.TrialDifferencePercent), ...
    double(LegacyStats.TrialUnifiedMaxTFWERP),trialMaxP, ...
    double(LegacyStats.CircularDifferencePercent),circularDifference, ...
    circularDifference-double(LegacyStats.CircularDifferencePercent), ...
    double(LegacyStats.CircularUnifiedMaxTFWERP),circMaxP, ...
    legacyBoth,mirrorBoth,status, ...
    'VariableNames',{'Band','ROI_A','ROI_B', ...
    'LockedTrialDifferencePercent','MirrorTrialDifferencePercent', ...
    'TrialDifferenceChangePP','LockedTrialUnifiedMaxTFWERP', ...
    'MirrorTrialUnifiedMaxTFWERP','LockedCircularDifferencePercent', ...
    'MirrorCircularDifferencePercent','CircularDifferenceChangePP', ...
    'LockedCircularUnifiedMaxTFWERP','MirrorCircularUnifiedMaxTFWERP', ...
    'LockedBothUnifiedMaxT05','MirrorBothUnifiedMaxT05','CandidateStatus'});

CandidateSummary = makeCandidateSummary(Comparison,bands);
StabilityQC = subjectEffectStability(Legacy,SubjectEffects);

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
runID = ['BURST_EDGE_FAMILY_MIRROR_PADDED_' stamp];
runDir = fullfile(outputRoot,runID);
mkdir(runDir);
InputQC = table(nSubjects,numel(rois),numel(bands),nEdges,nTests,2, ...
    fs,nSamples,double(p.Results.nPermutations),double(p.Results.seed), ...
    height(CE),height(CI),all(tf), ...
    'VariableNames',{'NSubjects','NROIs','NBands','NEdgesPerBand', ...
    'NTestsPerNull','ZThreshold','SamplingRateHz','WindowSamples', ...
    'NPermutations','Seed','NCorticalEvents','NCorticalInventoryRows', ...
    'LockedEdgeKeysIdentical'});

writetable(InputQC,fullfile(runDir,'INPUT_QC.csv'));
writetable(MirrorStats, ...
    fullfile(runDir,'MIRROR_PADDED_UNIFIED_630_EDGE_STATS.csv'));
writetable(SubjectEffects, ...
    fullfile(runDir,'MIRROR_PADDED_SUBJECT_LEVEL_630_EDGE_EFFECTS.csv'));
writetable(Comparison,fullfile(runDir,'EDGE_SET_COMPARISON.csv'));
writetable(CandidateSummary,fullfile(runDir,'CANDIDATE_SUMMARY.csv'));
writetable(StabilityQC,fullfile(runDir,'SUBJECT_EFFECT_STABILITY_QC.csv'));
writetable(MirrorStats(mirrorBoth,:), ...
    fullfile(runDir,'MIRROR_PADDED_EDGES_SURVIVING_BOTH.csv'));
writetable(Comparison(legacyBoth,:), ...
    fullfile(runDir,'LOCKED_CANDIDATES_UNDER_MIRROR_PADDING.csv'));
writetable(Comparison(legacyBoth | mirrorBoth,:), ...
    fullfile(runDir,'CANDIDATE_UNION_FOR_ORTHOGONALIZATION.csv'));

Results = struct;
Results.cfg = struct('runID',runID,'runDir',runDir, ...
    'mirrorBurstFile',mirrorBurstFile, ...
    'legacyUnifiedFile',legacyUnifiedFile,'created',stamp, ...
    'nPermutations',double(p.Results.nPermutations), ...
    'seed',double(p.Results.seed));
Results.InputQC = InputQC;
Results.MirrorStats = MirrorStats;
Results.SubjectEffects = SubjectEffects;
Results.Comparison = Comparison;
Results.CandidateSummary = CandidateSummary;
Results.StabilityQC = StabilityQC;
save(fullfile(runDir,[runID '_RESULTS.mat']),'Results','-v7.3');

fprintf('\nMirror-padded 630-edge family complete\n%s\n',runDir);
disp(CandidateSummary)
fprintf('\nLOCKED CANDIDATES UNDER MIRROR PADDING\n')
disp(Comparison(legacyBoth,:))
fprintf('\nCANDIDATE UNION FOR NEXT ORTHOGONALIZATION STEP\n')
disp(Comparison(legacyBoth | mirrorBoth,:))
end

function nSamples = inferWindowSamples(CI,CE)
if ismember('NSamples',CI.Properties.VariableNames)
    values = unique(double(CI.NSamples(isfinite(double(CI.NSamples)))));
    assert(numel(values)==1,'CorticalInventory contains inconsistent NSamples.');
    nSamples = values;
else
    nSamples = max(double(CE.OffsetSample));
end
end

function [trialIDs,masks] = makeMasks(CI,CE,sid,roi,band,condition,nSamples)
q = string(CI.Subject)==sid & string(CI.Signal)==roi & ...
    string(CI.Band)==band & string(CI.Condition)==condition & ...
    double(CI.ZThreshold)==2;
if ismember('TrialValid',CI.Properties.VariableNames)
    q = q & logical(CI.TrialValid);
end
trialIDs = unique(double(CI.OriginalTrial(q)),'stable');
assert(~isempty(trialIDs),'No valid trials for %s %s %s %s.', ...
    sid,roi,band,condition);
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
assert(sum(hits)==1, ...
    'Expected exactly one burst object with events/inventory.');
B = S.(names{find(hits,1)});
end

function U = loadUnifiedObject(filePath)
S = load(filePath);
names = fieldnames(S);
hits = false(numel(names),1);
for k = 1:numel(names)
    x = S.(names{k});
    hits(k) = isstruct(x) && isfield(x,'UnifiedStats');
end
assert(sum(hits)==1, ...
    'Expected one cropped-envelope object containing UnifiedStats.');
U = S.(names{find(hits,1)});
end

function key = edgeKey(band,a,b)
band = string(band);
a = string(a);
b = string(b);
lo = strings(size(a));
hi = strings(size(a));
for k = 1:numel(a)
    pair = sort([a(k) b(k)]);
    lo(k) = pair(1);
    hi(k) = pair(2);
end
key = band + "|" + lo + "|" + hi;
end

function Summary = makeCandidateSummary(C,bands)
Band = [bands;"ALL"];
n = numel(Band);
Locked = zeros(n,1);
Mirror = zeros(n,1);
Retained = zeros(n,1);
Lost = zeros(n,1);
Gained = zeros(n,1);
for k = 1:n
    if Band(k)=="ALL"
        q = true(height(C),1);
    else
        q = string(C.Band)==Band(k);
    end
    Locked(k) = sum(C.LockedBothUnifiedMaxT05(q));
    Mirror(k) = sum(C.MirrorBothUnifiedMaxT05(q));
    Retained(k) = sum(C.CandidateStatus(q)=="RETAINED_AFTER_MIRROR_PADDING");
    Lost(k) = sum(C.CandidateStatus(q)=="LOST_AFTER_MIRROR_PADDING");
    Gained(k) = sum(C.CandidateStatus(q)=="GAINED_AFTER_MIRROR_PADDING");
end
Summary = table(Band,Locked,Mirror,Retained,Lost,Gained, ...
    'VariableNames',{'Band','NLockedCandidates','NMirrorCandidates', ...
    'NRetained','NLost','NGained'});
end

function QC = subjectEffectStability(Legacy,MirrorEffects)
if ~isfield(Legacy,'SubjectEffects') || isempty(Legacy.SubjectEffects)
    QC = table("NOT_AVAILABLE",NaN,NaN,NaN,NaN,0, ...
        'VariableNames',{'Status','TrialPearsonR','CircularPearsonR', ...
        'TrialMeanAbsoluteChangePP','CircularMeanAbsoluteChangePP', ...
        'RowsCompared'});
    return
end
L = Legacy.SubjectEffects;
required = ["Subject" "Band" "ROI_A" "ROI_B" ...
    "TrialDifferencePercent" "CircularDifferencePercent"];
assert(all(ismember(required,string(L.Properties.VariableNames))), ...
    'Locked SubjectEffects lacks required columns.');
keyL = string(L.Subject)+"|"+edgeKey(L.Band,L.ROI_A,L.ROI_B);
keyM = string(MirrorEffects.Subject)+"|"+ ...
    edgeKey(MirrorEffects.Band,MirrorEffects.ROI_A,MirrorEffects.ROI_B);
[tf,loc] = ismember(keyM,keyL);
assert(all(tf) && numel(unique(loc))==height(MirrorEffects), ...
    'Locked and mirror-padded subject-effect keys do not match.');
lt = double(L.TrialDifferencePercent(loc));
lc = double(L.CircularDifferencePercent(loc));
mt = double(MirrorEffects.TrialDifferencePercent);
mc = double(MirrorEffects.CircularDifferencePercent);
QC = table("COMPLETE",pearsonManual(lt,mt),pearsonManual(lc,mc), ...
    mean(abs(mt-lt)),mean(abs(mc-lc)),height(MirrorEffects), ...
    'VariableNames',{'Status','TrialPearsonR','CircularPearsonR', ...
    'TrialMeanAbsoluteChangePP','CircularMeanAbsoluteChangePP', ...
    'RowsCompared'});
end

function r = pearsonManual(x,y)
x = double(x(:));
y = double(y(:));
q = isfinite(x) & isfinite(y);
x = x(q)-mean(x(q));
y = y(q)-mean(y(q));
denominator = sqrt(sum(x.^2)*sum(y.^2));
if denominator==0
    r = NaN;
else
    r = sum(x.*y)/denominator;
end
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
