function Results = burst_boundary_count_sensitivity(burstFile,varargin)
% Conservative sensitivity analysis for right-window-censored events.
%
% Recomputes the primary z=2 theta/alpha/beta count outcome after excluding
% all events whose detected component reaches the final sample of the
% 600-sample window. This does not redetect events and therefore does not
% replace the primary onset-observed count analysis. It tests whether the
% reported count contrasts depend on events touching the right boundary.

p = inputParser;
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
addParameter(p,'seed',20260921,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun,'Set confirmRun=true after checking the input path.');

burstFile = char(string(burstFile));
outputRoot = char(string(p.Results.outputRoot));
assert(isfile(burstFile),'Missing burst file: %s',burstFile);
if isempty(outputRoot), outputRoot = fileparts(burstFile); end

B = loadBurstObject(burstFile);
CE = B.CorticalEvents;
CI = B.CorticalInventory;

requiredEvent = ["Subject" "Signal" "Band" "Condition" "ZThreshold" ...
    "OriginalTrial" "PrimaryEligible" "RightCensored"];
requiredInventory = ["Subject" "Signal" "Band" "Condition" ...
    "ZThreshold" "OriginalTrial"];
assert(all(ismember(requiredEvent,string(CE.Properties.VariableNames))), ...
    'CorticalEvents lacks required columns.');
assert(all(ismember(requiredInventory,string(CI.Properties.VariableNames))), ...
    'CorticalInventory lacks required columns.');

subjects = sort(unique(string(CI.Subject)));
bands = ["theta";"alpha";"beta"];
conditions = ["EMG2";"EMG3"];
assert(numel(subjects)==73,'Expected 73 healthy participants.');

nRows = numel(subjects)*numel(bands)*numel(conditions);
R = cell(nRows,10); row = 0;
for s = 1:numel(subjects)
    sid = subjects(s);
    for b = 1:numel(bands)
        band = bands(b);
        for c = 1:numel(conditions)
            condition = conditions(c);
            iq = string(CI.Subject)==sid & string(CI.Band)==band & ...
                string(CI.Condition)==condition & double(CI.ZThreshold)==2;
            if ismember('TrialValid',CI.Properties.VariableNames)
                iq = iq & logical(CI.TrialValid);
            end
            nValidROITrials = sum(iq);
            assert(nValidROITrials>0,'No valid ROI-trials for %s %s %s.', ...
                sid,band,condition);

            eq = string(CE.Subject)==sid & string(CE.Band)==band & ...
                string(CE.Condition)==condition & double(CE.ZThreshold)==2 & ...
                logical(CE.PrimaryEligible);
            nPrimary = sum(eq);
            nRight = sum(eq & logical(CE.RightCensored));
            nComplete = sum(eq & ~logical(CE.RightCensored));
            assert(nPrimary==nRight+nComplete,'Boundary-event accounting failed.');

            % Reproduce the locked ROI-first aggregation exactly: calculate
            % the mean event count per valid trial separately for each ROI,
            % then take the median of the 21 ROI-specific values.
            roiNames = sort(unique(string(CI.Signal(iq))));
            assert(numel(roiNames)==21,'Expected 21 cortical ROIs.');
            primaryROI = nan(21,1);
            completeROI = nan(21,1);
            for r = 1:21
                roi = roiNames(r);
                roiInventory = iq & string(CI.Signal)==roi;
                nROITrials = sum(roiInventory);
                assert(nROITrials>0,'No valid trials for %s %s %s %s.', ...
                    sid,band,condition,roi);
                roiEvents = eq & string(CE.Signal)==roi;
                primaryROI(r) = sum(roiEvents)/nROITrials;
                completeROI(r) = sum(roiEvents & ...
                    ~logical(CE.RightCensored))/nROITrials;
            end
            primaryFeature = median(primaryROI,'omitnan');
            completeFeature = median(completeROI,'omitnan');

            row = row+1;
            R(row,:) = {sid,band,condition,nValidROITrials,nPrimary,nRight, ...
                nComplete,primaryFeature,completeFeature, ...
                100*nRight/max(nPrimary,1)};
        end
    end
end

ParticipantCondition = cell2table(R,'VariableNames',{ ...
    'Subject','Band','Condition','NValidROITrials','NPrimaryEvents', ...
    'NRightCensoredEvents','NCompleteEvents','PrimaryCountPerROITrial', ...
    'BoundaryExcludedCountPerROITrial','RightCensoredPercentOfPrimaryEvents'});

Dprimary = nan(numel(subjects),numel(bands));
Dcomplete = nan(size(Dprimary));
PrimaryMeanEMG2 = nan(numel(bands),1); PrimaryMeanEMG3 = PrimaryMeanEMG2;
CompleteMeanEMG2 = PrimaryMeanEMG2; CompleteMeanEMG3 = PrimaryMeanEMG2;
NRightEMG2 = zeros(numel(bands),1); NRightEMG3 = NRightEMG2;

for b = 1:numel(bands)
    band = bands(b);
    A = sortrows(ParticipantCondition( ...
        ParticipantCondition.Band==band & ...
        ParticipantCondition.Condition=="EMG2",:),'Subject');
    C = sortrows(ParticipantCondition( ...
        ParticipantCondition.Band==band & ...
        ParticipantCondition.Condition=="EMG3",:),'Subject');
    assert(isequal(A.Subject,C.Subject) && isequal(A.Subject,subjects), ...
        'Participant alignment failed for %s.',band);
    Dprimary(:,b) = A.PrimaryCountPerROITrial-C.PrimaryCountPerROITrial;
    Dcomplete(:,b) = A.BoundaryExcludedCountPerROITrial- ...
        C.BoundaryExcludedCountPerROITrial;
    PrimaryMeanEMG2(b) = mean(A.PrimaryCountPerROITrial);
    PrimaryMeanEMG3(b) = mean(C.PrimaryCountPerROITrial);
    CompleteMeanEMG2(b) = mean(A.BoundaryExcludedCountPerROITrial);
    CompleteMeanEMG3(b) = mean(C.BoundaryExcludedCountPerROITrial);
    NRightEMG2(b) = sum(A.NRightCensoredEvents);
    NRightEMG3(b) = sum(C.NRightCensoredEvents);
end

% Broad reproduction gate against the rounded primary values in Table 1.
expectedEMG2 = [0.0038;0.0676;0.0467];
expectedEMG3 = [0.0007;0.0504;0.0247];
reproductionPass = all(abs(PrimaryMeanEMG2-expectedEMG2)<5e-4) && ...
    all(abs(PrimaryMeanEMG3-expectedEMG3)<5e-4);
assert(reproductionPass, ...
    'Primary count means do not reproduce the locked Table 1 values.');

rng(double(p.Results.seed),'twister');
[tComplete,dzComplete,rawPComplete,maxPComplete] = ...
    maxT(Dcomplete,double(p.Results.nPermutations));
rng(double(p.Results.seed)+1,'twister');
ciComplete = bootstrapCI(Dcomplete,double(p.Results.nBootstrap));

PrimaryDifference = mean(Dprimary,1)';
BoundaryExcludedDifference = mean(Dcomplete,1)';
DifferenceChange = BoundaryExcludedDifference-PrimaryDifference;
RelativeEffectRetainedPercent = 100*BoundaryExcludedDifference./PrimaryDifference;

Summary = table(bands,repmat(73,3,1),PrimaryMeanEMG2,PrimaryMeanEMG3, ...
    PrimaryDifference,CompleteMeanEMG2,CompleteMeanEMG3, ...
    BoundaryExcludedDifference,ciComplete(:,1),ciComplete(:,2), ...
    tComplete,dzComplete,rawPComplete,maxPComplete, ...
    BoundaryExcludedDifference>0 & maxPComplete<0.05, ...
    DifferenceChange,RelativeEffectRetainedPercent,NRightEMG2,NRightEMG3, ...
    'VariableNames',{'Band','N','PrimaryMeanEMG2','PrimaryMeanEMG3', ...
    'PrimaryDifference','BoundaryExcludedMeanEMG2','BoundaryExcludedMeanEMG3', ...
    'BoundaryExcludedDifference','BootstrapCI95Lower','BootstrapCI95Upper', ...
    'T','CohenDz','RawSignFlipP','ThreeBandMaxTFWERP','PositiveAndMaxT05', ...
    'DifferenceChange','RelativeEffectRetainedPercent', ...
    'NRightCensoredEventsEMG2','NRightCensoredEventsEMG3'});

SubjectDifferences = table(subjects,Dprimary(:,1),Dprimary(:,2),Dprimary(:,3), ...
    Dcomplete(:,1),Dcomplete(:,2),Dcomplete(:,3), ...
    'VariableNames',{'Subject','PrimaryThetaDifference','PrimaryAlphaDifference', ...
    'PrimaryBetaDifference','BoundaryExcludedThetaDifference', ...
    'BoundaryExcludedAlphaDifference','BoundaryExcludedBetaDifference'});

stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
runID = ['BURST_BOUNDARY_COUNT_SENSITIVITY_' stamp];
runDir = fullfile(outputRoot,runID);
[ok,msg] = mkdir(runDir); assert(ok,'Cannot create %s: %s',runDir,msg);

InputQC = table(73,21,3,2,600,2,double(p.Results.nPermutations), ...
    double(p.Results.nBootstrap),double(p.Results.seed),reproductionPass, ...
    sum(NRightEMG2)+sum(NRightEMG3), ...
    'VariableNames',{'NSubjects','NROIs','NBands','NConditions', ...
    'WindowSamples','ZThreshold','NPermutations','NBootstrap','Seed', ...
    'PrimaryTableMeansReproduced','NRightCensoredStandardBandEvents'});

writetable(InputQC,fullfile(runDir,'INPUT_QC.csv'));
writetable(ParticipantCondition,fullfile(runDir,'PARTICIPANT_CONDITION_COUNTS.csv'));
writetable(SubjectDifferences,fullfile(runDir,'SUBJECT_DIFFERENCES.csv'));
writetable(Summary,fullfile(runDir,'RIGHT_BOUNDARY_COUNT_SENSITIVITY.csv'));

Results = struct;
Results.cfg = struct('runID',runID,'runDir',runDir,'burstFile',burstFile, ...
    'created',stamp,'nPermutations',double(p.Results.nPermutations), ...
    'nBootstrap',double(p.Results.nBootstrap),'seed',double(p.Results.seed), ...
    'countAggregation', ...
    'median across 21 ROI-specific mean event counts per valid trial', ...
    'interpretation',['Conservative sensitivity excluding detected events ' ...
    'that touch the final sample; no burst redetection was performed.']);
Results.InputQC = InputQC;
Results.ParticipantCondition = ParticipantCondition;
Results.SubjectDifferences = SubjectDifferences;
Results.Summary = Summary;
save(fullfile(runDir,[runID '_RESULTS.mat']),'Results','-v7.3');

fprintf('\nRight-boundary count sensitivity complete:\n%s\n',runDir);
disp(InputQC)
disp(Summary)
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

function ci = bootstrapCI(D,nBootstrap)
[n,m] = size(D); ci = nan(m,2);
for j = 1:m
    idx = randi(n,n,nBootstrap);
    sampled = reshape(D(idx(:),j),n,nBootstrap);
    bootMeans = mean(sampled,1);
    ci(j,:) = prctile(bootMeans,[2.5 97.5]);
end
end

function [obsT,dz,rawP,maxP] = maxT(D,nPerm)
[n,m] = size(D);
mu = mean(D,1); sd = std(D,0,1);
valid = isfinite(mu) & isfinite(sd) & sd>0;
obsT = zeros(m,1); dz = zeros(m,1);
obsT(valid) = (mu(valid)./(sd(valid)/sqrt(n)))';
dz(valid) = (mu(valid)./sd(valid))';
rawCount = zeros(1,m); maxCount = zeros(1,m);
sumSq = sum(D.^2,1); batchSize = 1000; done = 0;
while done<nPerm
    batch = min(batchSize,nPerm-done);
    signs = 2*(rand(batch,n)>0.5)-1;
    permMean = (signs*D)/n;
    numerator = max(sumSq-n*(permMean.^2),0);
    permSD = sqrt(numerator/(n-1));
    permT = zeros(batch,m);
    permT(:,valid) = permMean(:,valid)./(permSD(:,valid)/sqrt(n));
    rawCount = rawCount+sum(abs(permT)>=abs(obsT'),1);
    maximumT = max(abs(permT),[],2);
    maxCount = maxCount+sum(maximumT>=abs(obsT'),1);
    done = done+batch;
end
rawP = ((rawCount+1)/(nPerm+1))';
maxP = ((maxCount+1)/(nPerm+1))';
rawP(~valid) = 1; maxP(~valid) = 1;
end
