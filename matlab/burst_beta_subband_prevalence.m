function [SubbandStats,cfgStats] = burst_beta_subband_prevalence( ...
    burstFullFile,varargin)
% Prespecified participant-level inference for low-/high-beta primary burst
% prevalence. The primary z=2 family comprises low-beta EMG2-EMG3,
% high-beta EMG2-EMG3, and their difference-of-differences. Two-sided
% sign-flip tests use joint max-|T| FWER correction over the three tests.

p=inputParser;
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
addParameter(p,'seed',42008,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun, ...
    'Set confirmRun=true only after accepting the locked analysis family.');

burstFullFile=char(string(burstFullFile));
assert(isfile(burstFullFile),'Burst full-cohort file is missing.');
Z=load(burstFullFile,'BurstFull','cfgBurst');
B=Z.BurstFull;cfgBurst=Z.cfgBurst;
I=B.CorticalInventory;
requiredVars=["Subject" "Signal" "Band" "Condition" ...
    "ZThreshold" "AnyBurstPrimary"];
assert(all(ismember(requiredVars,string(I.Properties.VariableNames))), ...
    'CorticalInventory lacks required variables.');
assert(all(B.NestingQC.NViolationsInclusive==0)& ...
    all(B.NestingQC.NViolationsPrimary==0),'Cortical nesting QC failed.');
assert(sum(B.CorticalNormalization.ScaleFallback)==0, ...
    'Cortical normalization fallback detected.');

subjects=string(cfgBurst.subjects(:));
bands=["lowbeta";"highbeta"];
conditions=["EMG2";"EMG3"];
thresholds=[1.5 2 2.5];
assert(numel(subjects)==73&&numel(unique(subjects))==73, ...
    'Expected 73 unique healthy participants.');
assert(all(ismember(bands,string(cfgBurst.bandLabels(:)))), ...
    'Required beta subbands are absent.');

% One network-averaged prevalence per participant, band, condition and
% threshold. Each value is the mean across all valid ROI-trials in the
% corresponding cell; participants, not ROI-trials, are the inference unit.
X=nan(numel(subjects),numel(bands),numel(conditions),numel(thresholds));
rows=cell(numel(subjects)*numel(bands)*numel(conditions)*numel(thresholds),6);
r=0;
for s=1:numel(subjects)
    for b=1:numel(bands)
        for c=1:numel(conditions)
            for z=1:numel(thresholds)
                q=string(I.Subject)==subjects(s)& ...
                    string(I.Band)==bands(b)& ...
                    string(I.Condition)==conditions(c)& ...
                    I.ZThreshold==thresholds(z);
                assert(any(q),'Missing participant cell: %s %s %s z=%g', ...
                    subjects(s),bands(b),conditions(c),thresholds(z));
                value=mean(double(I.AnyBurstPrimary(q)));
                assert(isfinite(value)&&value>=0&&value<=1, ...
                    'Invalid prevalence value.');
                X(s,b,c,z)=value;
                r=r+1;
                rows(r,:)={subjects(s),bands(b),conditions(c),thresholds(z), ...
                    sum(q),value};
            end
        end
    end
end
ParticipantPrevalence=cell2table(rows,'VariableNames',{ ...
    'Subject','Band','Condition','ZThreshold','NValidROITrials', ...
    'PrimaryBurstPrevalence'});
ParticipantPrevalence.Subject=string(ParticipantPrevalence.Subject);
ParticipantPrevalence.Band=string(ParticipantPrevalence.Band);
ParticipantPrevalence.Condition=string(ParticipantPrevalence.Condition);

allParts=cell(numel(thresholds),1);
for z=1:numel(thresholds)
    lowDifference=100*(X(:,1,1,z)-X(:,1,2,z));
    highDifference=100*(X(:,2,1,z)-X(:,2,2,z));
    D=[lowDifference highDifference highDifference-lowDifference];
    contrast=["low beta: EMG2 - EMG3"; ...
        "high beta: EMG2 - EMG3"; ...
        "interaction: high-beta difference - low-beta difference"];
    allParts{z}=runFamily(D,contrast,thresholds(z), ...
        p.Results.nPermutations,p.Results.nBootstrap,p.Results.seed+z-1);
end
AllThresholdContrasts=vertcat(allParts{:});
PrimaryContrasts=AllThresholdContrasts(AllThresholdContrasts.ZThreshold==2,:);
ThresholdSensitivity=AllThresholdContrasts(AllThresholdContrasts.ZThreshold~=2,:);

outputRoot=char(string(p.Results.outputRoot));
if isempty(outputRoot),outputRoot=char(cfgBurst.runDir);end
timeTag=char(datetime('now','Format','yyyyMMdd_HHmmss'));
runID=['BURST_BETA_SUBBAND_PREVALENCE_' timeTag];
runDir=fullfile(outputRoot,runID);
assert(~isfolder(runDir),'Refusing overwrite: %s',runDir);
[ok,msg]=mkdir(runDir);assert(ok,'Cannot create %s: %s',runDir,msg);
cfgStats=struct('runID',runID,'runDir',string(runDir), ...
    'burstFullFile',string(burstFullFile),'sourceBurstRunID',cfgBurst.runID, ...
    'subjects',subjects,'bands',bands,'conditions',conditions, ...
    'primaryZThreshold',2,'sensitivityZThresholds',[1.5 2.5], ...
    'nPermutations',p.Results.nPermutations, ...
    'nBootstrap',p.Results.nBootstrap,'seed',p.Results.seed, ...
    'inferenceUnit','participant', ...
    'primaryFamily','three prevalence contrasts; joint max-absolute-T FWER', ...
    'analysisStatus','EXPLORATORY_PRESPECIFIED_BEFORE_FULLCOHORT_RESULTS');
SubbandStats=struct('ParticipantPrevalence',ParticipantPrevalence, ...
    'PrimaryContrasts',PrimaryContrasts, ...
    'ThresholdSensitivity',ThresholdSensitivity, ...
    'AllThresholdContrasts',AllThresholdContrasts);
save(fullfile(runDir,[runID '_RESULTS.mat']),'SubbandStats','cfgStats','-v7.3');
writetable(ParticipantPrevalence,fullfile(runDir,'PARTICIPANT_PREVALENCE.csv'));
writetable(PrimaryContrasts,fullfile(runDir,'PRIMARY_Z2_CONTRASTS.csv'));
writetable(ThresholdSensitivity,fullfile(runDir,'THRESHOLD_SENSITIVITY.csv'));
fprintf('\nBeta-subband prevalence analysis complete:\n%s\n',runDir);
fprintf('\nPRIMARY z=2 FAMILY\n');disp(PrimaryContrasts)
fprintf('\nTHRESHOLD SENSITIVITY\n');disp(ThresholdSensitivity)
end

function T=runFamily(D,contrast,zThreshold,nPerm,nBoot,seed)
assert(size(D,2)==3&&all(isfinite(D(:))), ...
    'Expected a complete N x 3 contrast matrix.');
n=size(D,1);assert(n==73,'Expected 73 participants.');
mu=mean(D,1);sd=std(D,0,1);
assert(all(sd>0),'At least one contrast has zero variance.');
tObserved=mu./(sd/sqrt(n));
cohenDz=mu./sd;
proportionPositive=mean(D>0,1);

rng(seed,'twister');
rawExceed=zeros(1,3);maxTExceed=zeros(1,3);
chunkSize=5000;completed=0;
while completed<nPerm
    m=min(chunkSize,nPerm-completed);
    signs=2*randi([0 1],n,m)-1;
    permT=zeros(m,3);
    for j=1:3
        permMean=(signs'*D(:,j))/n;
        variance=(sum(D(:,j).^2)-n*permMean.^2)/(n-1);
        variance=max(variance,0);
        permT(:,j)=permMean./(sqrt(variance)/sqrt(n));
    end
    maxAbsT=max(abs(permT),[],2);
    for j=1:3
        rawExceed(j)=rawExceed(j)+sum(abs(permT(:,j))>=abs(tObserved(j)));
        maxTExceed(j)=maxTExceed(j)+sum(maxAbsT>=abs(tObserved(j)));
    end
    completed=completed+m;
end
rawP=(rawExceed+1)/(nPerm+1);
maxTFWERP=(maxTExceed+1)/(nPerm+1);

rng(seed+1000,'twister');
bootIndex=randi(n,n,nBoot);
bootMean=zeros(nBoot,3);
for j=1:3
    values=D(:,j);
    bootMean(:,j)=mean(values(bootIndex),1)';
end
ci=prctile(bootMean,[2.5 97.5],1);

T=table(contrast(:),repmat(zThreshold,3,1),repmat(n,3,1), ...
    mu(:),ci(1,:)',ci(2,:)',cohenDz(:),tObserved(:),rawP(:), ...
    maxTFWERP(:),proportionPositive(:), ...
    'VariableNames',{'Contrast','ZThreshold','NSubjects', ...
    'MeanDifferencePP','BootstrapCI95Lower','BootstrapCI95Upper', ...
    'CohenDz','T','RawP','MaxTFWERP','ProportionSubjectsPositive'});
end
