function Results = burst_orthogonalized_edges(varargin)
% Pairwise bidirectional orthogonalization sensitivity for the union of
% cropped-envelope and mirror-padded co-burst candidates. Analytic signals
% are calculated after symmetric reflection padding. The eight candidate
% edges are corrected jointly for each null-derived contrast.

p = inputParser;
addParameter(p,'sourceRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'candidateFile','',@(x)ischar(x)||isstring(x));
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'mirrorPadSamples',599,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'zThreshold',2,@(x)isnumeric(x)&&isscalar(x)&&x==2);
addParameter(p,'parentThreshold',1.5,@(x)isnumeric(x)&&isscalar(x)&&x==1.5);
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
addParameter(p,'randomSeed',20260924,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun, ...
    'Set confirmRun=true only after checking all input paths.');

sourceRoot = char(string(p.Results.sourceRoot));
candidateFile = char(string(p.Results.candidateFile));
outputRoot = char(string(p.Results.outputRoot));
assert(isfolder(sourceRoot),'Source root missing: %s',sourceRoot);
assert(isfile(candidateFile),'Candidate-union file missing: %s',candidateFile);
assert(isfolder(outputRoot),'Output root missing: %s',outputRoot);

Candidates = readtable(candidateFile,'TextType','string');
required = ["Band" "ROI_A" "ROI_B" "LockedBothUnifiedMaxT05" ...
    "MirrorBothUnifiedMaxT05" "CandidateStatus"];
assert(all(ismember(required,string(Candidates.Properties.VariableNames))), ...
    'Candidate-union table lacks required columns.');
Candidates.Band = string(Candidates.Band);
Candidates.ROI_A = string(Candidates.ROI_A);
Candidates.ROI_B = string(Candidates.ROI_B);
assert(height(Candidates)==8, ...
    'Expected the eight-edge locked/mirror candidate union; found %d.', ...
    height(Candidates));
assert(sum(Candidates.Band=="alpha")==2 && ...
    sum(Candidates.Band=="beta")==6 && ...
    ~any(Candidates.Band=="theta"), ...
    'Expected two alpha and six beta candidate edges.');
candidateKeys = edgeKey(Candidates.Band,Candidates.ROI_A,Candidates.ROI_B);
assert(numel(unique(candidateKeys))==height(Candidates), ...
    'Candidate-union table contains duplicate edges.');

conditions = ["EMG2";"EMG3"];
bands = ["alpha";"beta"];
bandHz = [8 13;13 30];
nSamples = 600;
fs = 600;
nSubjects = 73;
nEdges = height(Candidates);

allFiles = dir(fullfile(sourceRoot,'**','*_ROI_PCA_PILOT.mat'));
allPaths = strings(numel(allFiles),1);
for k = 1:numel(allFiles)
    allPaths(k) = string(fullfile(allFiles(k).folder,allFiles(k).name));
end
selected = false(numel(allFiles),1);
for b = 1:numel(bands)
    selected = selected | contains(upper(allPaths),"_"+upper(bands(b))+"_");
end
allFiles = allFiles(selected);
allPaths = allPaths(selected);
assert(numel(allFiles)==146, ...
    'Expected 73 alpha and 73 beta ROI-PCA files; found %d.',numel(allFiles));

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
runID = ['BURST_ORTHOGONALIZED_EDGES_' ...
    'CANDIDATE_UNION_R01_' stamp];
runDir = fullfile(outputRoot,runID);
assert(~isfolder(runDir),'Refusing overwrite: %s',runDir);
[ok,msg] = mkdir(runDir);
assert(ok,'Cannot create output directory: %s',msg);

subjectRows = cell(nSubjects*2*nEdges,26);
normalizationRows = cell(nSubjects*nEdges*4,9);
row = 0;
normRow = 0;
timer = tic;
subjectSet = strings(0,1);

for b = 1:numel(bands)
    band = bands(b);
    bandCandidates = Candidates(Candidates.Band==band,:);
    minimumDurationSamples = ceil(fs*2/mean(bandHz(b,:)));
    bandFileMask = contains(upper(allPaths),"_"+upper(band)+"_");
    bandFiles = allFiles(bandFileMask);
    assert(numel(bandFiles)==73,'Expected 73 %s ROI-PCA files.',band);
    fprintf('\nOrthogonalization: %s | %d candidate edges\n', ...
        upper(band),height(bandCandidates));

    bandSubjects = strings(numel(bandFiles),1);
    for s = 1:numel(bandFiles)
        file = fullfile(bandFiles(s).folder,bandFiles(s).name);
        X = load(file,'ROIPilot');
        assert(isfield(X,'ROIPilot'),'ROIPilot missing: %s',file);
        P = X.ROIPilot;
        sid = string(P.subject);
        bandSubjects(s) = sid;
        assert(string(P.bandLabel)==band,'Incorrect band file selected.');
        labels = string(P.roiLabels(:));
        assert(numel(labels)==21,'Expected 21 ROIs.');
        validateTrialKeys(P);

        for e = 1:height(bandCandidates)
            roiA = bandCandidates.ROI_A(e);
            roiB = bandCandidates.ROI_B(e);
            ia = find(labels==roiA);
            ib = find(labels==roiB);
            assert(isscalar(ia) && isscalar(ib), ...
                'Candidate ROI pair missing for %s: %s - %s.',sid,roiA,roiB);

            envelopes = cell(2,4);
            for ci = 1:2
                nTrial = size(P.roiData{ci},2);
                for q = 1:4
                    envelopes{ci,q} = nan(nTrial,nSamples);
                end
                for t = 1:nTrial
                    a = squeeze(double(P.roiData{ci}(ia,t,:)))';
                    bSignal = squeeze(double(P.roiData{ci}(ib,t,:)))';
                    assert(numel(a)==nSamples && numel(bSignal)==nSamples, ...
                        'Unexpected source time-series length.');
                    [za,zb] = mirrorPaddedAnalyticPair( ...
                        a,bSignal,double(p.Results.mirrorPadSamples));
                    ampA = abs(za);
                    ampB = abs(zb);
                    denomA = max(ampA,eps(max(1,max(ampA))));
                    denomB = max(ampB,eps(max(1,max(ampB))));
                    envelopes{ci,1}(t,:) = ampA;
                    envelopes{ci,2}(t,:) = ampB;
                    envelopes{ci,3}(t,:) = ...
                        abs(imag(za.*conj(zb)./denomB));
                    envelopes{ci,4}(t,:) = ...
                        abs(imag(zb.*conj(za)./denomA));
                end
            end

            centers = nan(1,4);
            scales = nan(1,4);
            fallback = false(1,4);
            definitions = ["A_original","B_original", ...
                "A_orthogonal_to_B","B_orthogonal_to_A"];
            for q = 1:4
                pooled = [envelopes{1,q}(:);envelopes{2,q}(:)];
                centers(q) = median(pooled,'omitnan');
                scales(q) = 1.4826*median(abs(pooled-centers(q)),'omitnan');
                if ~isfinite(scales(q)) || ...
                        scales(q)<=eps(max(1,abs(centers(q))))
                    scales(q) = std(pooled,0,'omitnan');
                    fallback(q) = true;
                end
                assert(isfinite(centers(q)) && isfinite(scales(q)) && ...
                    scales(q)>0,'Invalid normalization.');
                normRow = normRow+1;
                normalizationRows(normRow,:) = {sid,band,roiA,roiB, ...
                    definitions(q),centers(q),scales(q),numel(pooled), ...
                    fallback(q)};
            end

            for ci = 1:2
                masks = cell(1,4);
                for q = 1:4
                    z = (envelopes{ci,q}-centers(q))./scales(q);
                    masks{q} = primaryBurstMasks(z, ...
                        p.Results.parentThreshold,p.Results.zThreshold, ...
                        minimumDurationSamples);
                end
                D1 = edgeMetrics(masks{3},masks{2});
                D2 = edgeMetrics(masks{1},masks{4});
                symmetric = mean([D1;D2],1);
                row = row+1;
                subjectRows(row,:) = {sid,band,conditions(ci),roiA,roiB, ...
                    size(masks{1},1),p.Results.zThreshold, ...
                    minimumDurationSamples,D1(1),D1(2),D1(3),D1(4),D1(5), ...
                    D2(1),D2(2),D2(3),D2(4),D2(5), ...
                    symmetric(1),symmetric(2),symmetric(3), ...
                    symmetric(4),symmetric(5), ...
                    mean(masks{3},'all')*100, ...
                    mean(masks{4},'all')*100,file};
            end
        end
        fprintf('%s %d/73 | %s | elapsed %s\n',upper(band),s,sid, ...
            char(duration(0,0,toc(timer),'Format','hh:mm:ss')));
    end
    assert(numel(unique(bandSubjects))==73, ...
        'Duplicate participant files detected for %s.',band);
    if isempty(subjectSet)
        subjectSet = sort(bandSubjects);
    else
        assert(isequal(subjectSet,sort(bandSubjects)), ...
            'Alpha and beta participant sets differ.');
    end
end

assert(row==size(subjectRows,1) && normRow==size(normalizationRows,1), ...
    'Unexpected output row count.');
SubjectConditionEdge = cell2table(subjectRows,'VariableNames',{ ...
    'Subject','Band','Condition','ROI_A','ROI_B','NTrials','ZThreshold', ...
    'MinimumDurationSamples','D1ObservedOverlapPercent', ...
    'D1TrialNullPercent','D1TrialSpecificExcessPercent', ...
    'D1CircularNullPercent','D1TimeLockedExcessPercent', ...
    'D2ObservedOverlapPercent','D2TrialNullPercent', ...
    'D2TrialSpecificExcessPercent','D2CircularNullPercent', ...
    'D2TimeLockedExcessPercent','SymmetricObservedOverlapPercent', ...
    'SymmetricTrialNullPercent','SymmetricTrialSpecificExcessPercent', ...
    'SymmetricCircularNullPercent','SymmetricTimeLockedExcessPercent', ...
    'AOrthogonalizedOccupancyPercent','BOrthogonalizedOccupancyPercent', ...
    'SourceROIFile'});
Normalization = cell2table(normalizationRows,'VariableNames',{ ...
    'Subject','Band','ROI_A','ROI_B','SignalDefinition','NormCenter', ...
    'NormScale','NValues','ScaleFallback'});

[SubjectContrasts,EdgeStats] = inferContrasts( ...
    SubjectConditionEdge,Candidates,double(p.Results.nPermutations), ...
    double(p.Results.nBootstrap),double(p.Results.randomSeed));

% Carry the mirror-padded selection provenance into the final edge table.
keyStats = edgeKey(EdgeStats.Band,EdgeStats.ROI_A,EdgeStats.ROI_B);
[tf,loc] = ismember(keyStats,candidateKeys);
assert(all(tf) && numel(unique(loc))==nEdges, ...
    'Candidate provenance alignment failed.');
EdgeStats.LockedBothUnifiedMaxT05 = ...
    asLogical(Candidates.LockedBothUnifiedMaxT05(loc));
EdgeStats.MirrorBothUnifiedMaxT05 = ...
    asLogical(Candidates.MirrorBothUnifiedMaxT05(loc));
EdgeStats.MirrorPaddedCandidateStatus = string(Candidates.CandidateStatus(loc));

InputQC = table(nSubjects,21,2,nEdges,2,6, ...
    double(p.Results.zThreshold),double(p.Results.parentThreshold), ...
    double(p.Results.mirrorPadSamples),fs, ...
    double(p.Results.nPermutations),double(p.Results.nBootstrap), ...
    sum(Normalization.ScaleFallback),string(version), ...
    'VariableNames',{'NSubjects','NROIs','NBands','NCandidateEdges', ...
    'NAlphaCandidates','NBetaCandidates','ZThreshold','ParentThreshold', ...
    'MirrorPadSamples','SamplingRateHz','NPermutations','NBootstrap', ...
    'NNormalizationFallbacks','MATLABVersion'});
Decision = makeDecision(EdgeStats);
cfg = struct('runID',runID,'runDir',runDir,'sourceRoot',sourceRoot, ...
    'candidateFile',candidateFile, ...
    'method','bidirectional pairwise orthogonalization after mirror-padded analytic-signal calculation', ...
    'direction1','burst(A orthogonal to B) versus burst(original B)', ...
    'direction2','burst(original A) versus burst(B orthogonal to A)', ...
    'symmetricStatistic','arithmetic mean of the two directional overlap metrics', ...
    'normalization','joint EMG2+EMG3 median and 1.4826*MAD separately for each transformed signal', ...
    'multiplicity','separate max-|T| FWER across the eight-edge locked/mirror candidate union for each null-derived contrast', ...
    'interpretation','bounded leakage and edge-boundary sensitivity; candidate union is not an independent confirmatory family');
Results = struct('InputQC',InputQC,'Candidates',Candidates, ...
    'Normalization',Normalization, ...
    'SubjectConditionEdge',SubjectConditionEdge, ...
    'SubjectContrasts',SubjectContrasts,'EdgeStats',EdgeStats, ...
    'Decision',Decision,'cfg',cfg);

save(fullfile(runDir,[runID '_RESULTS.mat']),'Results','-v7.3');
writetable(InputQC,fullfile(runDir,'INPUT_QC.csv'));
writetable(Candidates,fullfile(runDir,'CANDIDATE_UNION_INPUT.csv'));
writetable(Normalization,fullfile(runDir,'NORMALIZATION.csv'));
writetable(SubjectConditionEdge, ...
    fullfile(runDir,'SUBJECT_CONDITION_EDGE.csv'));
writetable(SubjectContrasts,fullfile(runDir,'SUBJECT_CONTRASTS.csv'));
writetable(EdgeStats,fullfile(runDir,'EDGE_STATS.csv'));
writetable(EdgeStats(EdgeStats.BothNullsUnionMaxT05,:), ...
    fullfile(runDir,'EDGES_SURVIVING_ORTHOGONALIZATION.csv'));
writetable(Decision,fullfile(runDir,'DECISION.csv'));

fprintf('\nMirror-padded orthogonalized edge analysis complete\n%s\n',runDir);
disp(InputQC)
disp(EdgeStats)
disp(Decision)
end

function validateTrialKeys(P)
C = string(P.trialInventory.Condition);
for ci = 1:2
    q = find(C==string(P.conditions(ci)));
    assert(numel(q)==size(P.roiData{ci},2),'Trial inventory mismatch.');
    keys = double(P.trialInventory.OriginalTrial(q));
    assert(numel(unique(keys))==numel(keys), ...
        'Duplicate original trial keys.');
end
end

function [za,zb] = mirrorPaddedAnalyticPair(a,b,padRequested)
a = double(a(:)');
b = double(b(:)');
assert(isequal(size(a),size(b)) && all(isfinite(a)) && all(isfinite(b)), ...
    'Analytic-signal pair must be equal-sized and finite.');
n = numel(a);
padN = min(round(padRequested),n-1);
aPadded = [fliplr(a(1:padN)) a fliplr(a(end-padN+1:end))];
bPadded = [fliplr(b(1:padN)) b fliplr(b(end-padN+1:end))];
zaPadded = analyticSignalFFT(aPadded);
zbPadded = analyticSignalFFT(bPadded);
keep = padN+(1:n);
za = zaPadded(keep);
zb = zbPadded(keep);
end

function z = analyticSignalFFT(x)
x = double(x(:)');
n = numel(x);
multiplier = zeros(1,n);
if rem(n,2)==0
    multiplier(1) = 1;
    multiplier(n/2+1) = 1;
    multiplier(2:n/2) = 2;
else
    multiplier(1) = 1;
    multiplier(2:(n+1)/2) = 2;
end
z = ifft(fft(x).*multiplier);
end

function masks = primaryBurstMasks(z,parentThreshold,zThreshold,minSamples)
nTrial = size(z,1);
nSamples = size(z,2);
masks = false(nTrial,nSamples);
for t = 1:nTrial
    parents = findSegments(z(t,:)>=parentThreshold);
    parentID = zeros(1,nSamples);
    parentLeft = false(size(parents,1),1);
    for p = 1:size(parents,1)
        parentID(parents(p,1):parents(p,2)) = p;
        parentLeft(p) = parents(p,1)==1;
    end
    segments = findSegments(z(t,:)>=zThreshold);
    segments = segments((segments(:,2)-segments(:,1)+1)>=minSamples,:);
    for e = 1:size(segments,1)
        first = segments(e,1);
        last = segments(e,2);
        pid = parentID(first);
        assert(pid>0,'z=2 event lacks a z=1.5 parent.');
        if ~parentLeft(pid)
            masks(t,first:last) = true;
        end
    end
end
end

function M = edgeMetrics(A,B)
nTrial = size(A,1);
nSamples = size(A,2);
assert(isequal(size(A),size(B)) && nTrial>1,'Invalid paired masks.');
overlapByTrial = mean(A&B,2);
observed = mean(overlapByTrial);
pairOverlap = (double(A)*double(B)')/nSamples;
trialNull = (sum(pairOverlap,'all')-sum(diag(pairOverlap))) / ...
    (nTrial*(nTrial-1));
occupancyA = mean(A,2);
occupancyB = mean(B,2);
allShiftMean = occupancyA.*occupancyB;
nonzeroShiftMean = ...
    (nSamples*allShiftMean-overlapByTrial)/(nSamples-1);
circularNull = mean(nonzeroShiftMean);
M = 100*[observed,trialNull,observed-trialNull, ...
    circularNull,observed-circularNull];
end

function [S,E] = inferContrasts(T,Candidates,nPerm,nBoot,seed)
subjects = sort(unique(string(T.Subject)));
assert(numel(subjects)==73,'Expected 73 participants.');
nEdge = height(Candidates);
DT = nan(73,nEdge);
DC = nan(73,nEdge);
D1T = nan(73,nEdge);
D2T = nan(73,nEdge);
D1C = nan(73,nEdge);
D2C = nan(73,nEdge);
rows = cell(73*nEdge,11);
row = 0;
for e = 1:nEdge
    band = Candidates.Band(e);
    roiA = Candidates.ROI_A(e);
    roiB = Candidates.ROI_B(e);
    q = string(T.Band)==band & string(T.ROI_A)==roiA & ...
        string(T.ROI_B)==roiB;
    A = sortrows(T(q & string(T.Condition)=="EMG2",:),'Subject');
    B = sortrows(T(q & string(T.Condition)=="EMG3",:),'Subject');
    assert(height(A)==73 && height(B)==73 && ...
        isequal(string(A.Subject),subjects) && ...
        isequal(string(B.Subject),subjects),'Condition pairing failed.');
    DT(:,e) = A.SymmetricTrialSpecificExcessPercent- ...
        B.SymmetricTrialSpecificExcessPercent;
    DC(:,e) = A.SymmetricTimeLockedExcessPercent- ...
        B.SymmetricTimeLockedExcessPercent;
    D1T(:,e) = A.D1TrialSpecificExcessPercent-B.D1TrialSpecificExcessPercent;
    D2T(:,e) = A.D2TrialSpecificExcessPercent-B.D2TrialSpecificExcessPercent;
    D1C(:,e) = A.D1TimeLockedExcessPercent-B.D1TimeLockedExcessPercent;
    D2C(:,e) = A.D2TimeLockedExcessPercent-B.D2TimeLockedExcessPercent;
    for s = 1:73
        row = row+1;
        rows(row,:) = {subjects(s),band,roiA,roiB,DT(s,e),DC(s,e), ...
            D1T(s,e),D2T(s,e),D1C(s,e),D2C(s,e), ...
            sign(D1T(s,e))==sign(D2T(s,e))};
    end
end
S = cell2table(rows,'VariableNames',{'Subject','Band','ROI_A','ROI_B', ...
    'SymmetricTrialDifferencePercent','SymmetricCircularDifferencePercent', ...
    'D1TrialDifferencePercent','D2TrialDifferencePercent', ...
    'D1CircularDifferencePercent','D2CircularDifferencePercent', ...
    'TrialDirectionalSignAgreement'});

rng(seed,'twister');
[trialResult,trialCI] = maxTInference(DT,nPerm,nBoot);
rng(seed,'twister');
[circularResult,circularCI] = maxTInference(DC,nPerm,nBoot);
E = table(Candidates.Band,Candidates.ROI_A,Candidates.ROI_B, ...
    repmat(73,nEdge,1),mean(DT,1)',trialCI(:,1),trialCI(:,2), ...
    trialResult.T,trialResult.Dz,trialResult.RawP,trialResult.MaxTP, ...
    mean(DC,1)',circularCI(:,1),circularCI(:,2),circularResult.T, ...
    circularResult.Dz,circularResult.RawP,circularResult.MaxTP, ...
    mean(D1T,1)',mean(D2T,1)',mean(D1C,1)',mean(D2C,1)', ...
    mean(sign(D1T)==sign(D2T),1)', ...
    'VariableNames',{'Band','ROI_A','ROI_B','NSubjects', ...
    'TrialDifferencePercent','TrialCI95Lower','TrialCI95Upper', ...
    'TrialT','TrialCohenDz','TrialRawP','TrialUnionMaxTFWERP', ...
    'CircularDifferencePercent','CircularCI95Lower','CircularCI95Upper', ...
    'CircularT','CircularCohenDz','CircularRawP', ...
    'CircularUnionMaxTFWERP','D1TrialDifferencePercent', ...
    'D2TrialDifferencePercent','D1CircularDifferencePercent', ...
    'D2CircularDifferencePercent', ...
    'ProportionDirectionalTrialSignAgreement'});
E.TrialUnionMaxT05 = E.TrialDifferencePercent>0 & ...
    E.TrialUnionMaxTFWERP<0.05;
E.CircularUnionMaxT05 = E.CircularDifferencePercent>0 & ...
    E.CircularUnionMaxTFWERP<0.05;
E.BothNullsUnionMaxT05 = E.TrialUnionMaxT05 & E.CircularUnionMaxT05;
end

function [O,CI] = maxTInference(D,nPerm,nBoot)
n = size(D,1);
k = size(D,2);
mu = mean(D,1);
sd = std(D,0,1);
tObs = zeros(1,k);
valid = sd>0 & isfinite(sd);
tObs(valid) = mu(valid)./(sd(valid)/sqrt(n));
ss = sum(D.^2,1);
raw = zeros(1,k);
adjusted = zeros(1,k);
done = 0;
while done<nPerm
    m = min(2000,nPerm-done);
    signs = 2*(rand(m,n)>0.5)-1;
    nullMean = (signs*D)/n;
    nullVariance = (repmat(ss,m,1)-n*nullMean.^2)/(n-1);
    nullT = zeros(size(nullMean));
    good = nullVariance>0;
    nullT(good) = nullMean(good)./sqrt(nullVariance(good)/n);
    maximumT = max(abs(nullT),[],2);
    raw = raw+sum(abs(nullT)>=abs(tObs),1);
    adjusted = adjusted+sum(maximumT>=abs(tObs),1);
    done = done+m;
end
O = struct('T',tObs','Dz',(mu./sd)', ...
    'RawP',((raw+1)/(nPerm+1))','MaxTP',((adjusted+1)/(nPerm+1))');
O.Dz(~isfinite(O.Dz)) = 0;
CI = nan(k,2);
for j = 1:k
    CI(j,:) = bootstrapMeanCI(D(:,j),nBoot);
end
end

function ci = bootstrapMeanCI(x,nBoot)
n = numel(x);
boot = zeros(nBoot,1);
done = 0;
while done<nBoot
    m = min(2000,nBoot-done);
    idx = randi(n,n,m);
    sample = reshape(x(idx),n,m);
    boot(done+(1:m)) = mean(sample,1)';
    done = done+m;
end
ci = prctile(boot,[2.5 97.5]);
end

function D = makeDecision(E)
nTrial = sum(E.TrialUnionMaxT05);
nCircular = sum(E.CircularUnionMaxT05);
nBoth = sum(E.BothNullsUnionMaxT05);
status = "NO_EDGE_SURVIVES_ORTHOGONALIZED_BOTH_NULLS";
if nBoth==height(E)
    status = "ALL_UNION_EDGES_SURVIVE_ORTHOGONALIZED_BOTH_NULLS";
elseif nBoth>0
    status = "PARTIAL_EDGE_SURVIVAL_AFTER_ORTHOGONALIZATION";
end
D = table(status,nTrial,nCircular,nBoth,height(E), ...
    "Interpret as a conservative bounded sensitivity. Pairwise orthogonalization can remove genuine zero-lag neural coupling, and the locked/mirror candidate union is not an independent confirmatory family.", ...
    'VariableNames',{'Status','NTrialNullFWER05','NCircularNullFWER05', ...
    'NBothNullsFWER05','NUnionEdges','InterpretationGuard'});
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
key = band+"|"+lo+"|"+hi;
end

function q = asLogical(x)
if islogical(x)
    q = x;
elseif isnumeric(x)
    q = x~=0;
else
    value = lower(strtrim(string(x)));
    assert(all(ismember(value,["true" "false" "1" "0"])), ...
        'Cannot convert candidate flags to logical values.');
    q = value=="true" | value=="1";
end
q = logical(q(:));
end

function seg = findSegments(mask)
mask = logical(mask(:)');
d = diff([false mask false]);
seg = [find(d==1)' find(d==-1)'-1];
end
