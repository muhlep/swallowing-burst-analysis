function [DurationMatched,cfgMatched] = ...
    burst_beta_subband_duration_matched(burstFullFile,varargin)
% Sensitivity analysis imposing the low-beta absolute minimum duration
% (73 samples = 121.67 ms at 600 Hz) on both beta subbands. Source signals,
% normalization, thresholds and trial eligibility remain unchanged.

p=inputParser;
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'matchedMinimumSamples',73,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
addParameter(p,'minimumPairedN',20,@(x)isnumeric(x)&&isscalar(x)&&x>=10);
addParameter(p,'seed',42011,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun,'Set confirmRun=true after checking the input.');

burstFullFile=char(string(burstFullFile));assert(isfile(burstFullFile));
Z=load(burstFullFile,'BurstFull','cfgBurst');B=Z.BurstFull;cfgBurst=Z.cfgBurst;
I=B.CorticalInventory;E=B.CorticalEvents;
subjects=sort(unique(string(I.Subject)));bands=["lowbeta";"highbeta"];
conditions=["EMG2";"EMG3"];thresholds=[1.5 2 2.5];
matchedSamples=double(p.Results.matchedMinimumSamples);
assert(numel(subjects)==73&&matchedSamples==73, ...
    'This locked sensitivity requires 73 subjects and 73 samples.');
assert(all(B.NestingQC.NViolationsPrimary==0)&& ...
    sum(B.CorticalNormalization.ScaleFallback)==0,'Locked QC gate failed.');
requiredEventFields={'Subject','Signal','Band','Condition','OriginalTrial', ...
    'ZThreshold','DurationSamples','DurationSec','PeakZ','AUCAboveThreshold', ...
    'LeftCensored','RightCensored','PrimaryEligible'};
assert(all(ismember(requiredEventFields,E.Properties.VariableNames)), ...
    'Event fields missing.');

% Retain only events meeting the common absolute duration. Reconstruct the
% primary event count for every inventory key without changing any trial.
Em=E(double(E.DurationSamples)>=matchedSamples,:);
Ep=Em(logical(Em.PrimaryEligible),:);
eventKeys=makeKeys(Ep);[uniqueEventKeys,~,eventGroup]=unique(eventKeys);
eventCounts=accumarray(eventGroup,1);
inventoryKeys=makeKeys(I);[present,location]=ismember(inventoryKeys,uniqueEventKeys);
matchedCounts=zeros(height(I),1);matchedCounts(present)=eventCounts(location(present));
I.MatchedNBurstsPrimary=matchedCounts;
I.MatchedAnyBurstPrimary=matchedCounts>0;
lowRows=string(I.Band)=="lowbeta";
assert(all(I.MatchedNBurstsPrimary(lowRows)==I.NBurstsPrimary(lowRows)), ...
    'Low-beta reconstruction does not reproduce the original count.');
assert(all(I.MatchedNBurstsPrimary<=I.NBurstsPrimary), ...
    'Matched count exceeds original count.');
MatchQC=groupsummary(I,{'Band','Condition','ZThreshold'},{'mean','sum'}, ...
    {'AnyBurstPrimary','MatchedAnyBurstPrimary','NBurstsPrimary', ...
    'MatchedNBurstsPrimary'});

% Participant-level network prevalence at every threshold.
X=nan(73,2,2,3);
for s=1:73
    for b=1:2
        for c=1:2
            for z=1:3
                q=string(I.Subject)==subjects(s)&string(I.Band)==bands(b)& ...
                    string(I.Condition)==conditions(c)& ...
                    I.ZThreshold==thresholds(z)&I.TrialValid;
                assert(any(q),'Missing inventory cell.');
                X(s,b,c,z)=mean(double(I.MatchedAnyBurstPrimary(q)));
            end
        end
    end
end
prevalenceParts=cell(3,1);
for z=1:3
    low=100*(X(:,1,1,z)-X(:,1,2,z));
    high=100*(X(:,2,1,z)-X(:,2,2,z));
    D=[low high high-low];
    labels=["low beta: EMG2 - EMG3";"high beta: EMG2 - EMG3"; ...
        "interaction: high-beta difference - low-beta difference"];
    prevalenceParts{z}=runFamily(D,labels,repmat("percentage points",3,1), ...
        thresholds(z),p.Results.nPermutations,p.Results.nBootstrap, ...
        p.Results.minimumPairedN,p.Results.seed+z-1);
end
PrevalenceStats=vertcat(prevalenceParts{:});

% ROI-first feature aggregation at z=2 using duration-matched events.
metrics=["BurstCountPerROITrial";"DurationSec";"PeakZ";"AUCAboveThreshold"];
featureRows=cell(73*2*2*4,10);row=0;
for s=1:73
    sid=subjects(s);
    for b=1:2
        for c=1:2
            iq=string(I.Subject)==sid&string(I.Band)==bands(b)& ...
                string(I.Condition)==conditions(c)&I.ZThreshold==2&I.TrialValid;
            inv=I(iq,:);roiNames=unique(string(inv.Signal));
            assert(numel(roiNames)==21,'Expected 21 ROIs.');
            roiCount=nan(21,1);
            for r=1:21
                rq=string(inv.Signal)==roiNames(r);
                roiCount(r)=mean(double(inv.MatchedNBurstsPrimary(rq)));
            end
            row=row+1;featureRows(row,:)={sid,conditions(c),bands(b),metrics(1), ...
                median(roiCount,'omitnan'),21,height(inv), ...
                sum(inv.MatchedNBurstsPrimary),false,"ALL_VALID_ROI_TRIALS"};
            eq=string(Em.Subject)==sid&string(Em.Band)==bands(b)& ...
                string(Em.Condition)==conditions(c)&Em.ZThreshold==2& ...
                ~Em.LeftCensored&~Em.RightCensored;
            ev=Em(eq,:);fields={'DurationSec','PeakZ','AUCAboveThreshold'};
            for m=1:3
                roiValue=nan(21,1);nROIWithEvents=0;
                for r=1:21
                    rq=string(ev.Signal)==roiNames(r);
                    if any(rq)
                        nROIWithEvents=nROIWithEvents+1;
                        roiValue(r)=median(double(ev.(fields{m})(rq)),'omitnan');
                    end
                end
                row=row+1;featureRows(row,:)={sid,conditions(c),bands(b), ...
                    metrics(m+1),median(roiValue,'omitnan'),nROIWithEvents, ...
                    height(inv),height(ev),true,"COMPLETE_EVENTS_ROI_FIRST"};
            end
        end
    end
end
SubjectFeatures=cell2table(featureRows,'VariableNames',{'Subject','Condition', ...
    'Band','Metric','Value','NROIsContributing','NValidROITrials', ...
    'NObservations','ConditionalOnBurst','AggregationRule'});

D=nan(73,12);Metric=strings(12,1);Contrast=strings(12,1);
Unit=strings(12,1);endpoint=0;
units=["bursts per ROI-trial";"seconds";"z";"z-seconds"];
labels=["low beta: EMG2 - EMG3";"high beta: EMG2 - EMG3"; ...
    "interaction: high-beta difference - low-beta difference"];
for m=1:4
    V=nan(73,2,2);
    for b=1:2
        for c=1:2
            T=sortrows(SubjectFeatures(string(SubjectFeatures.Band)==bands(b)& ...
                string(SubjectFeatures.Condition)==conditions(c)& ...
                string(SubjectFeatures.Metric)==metrics(m),:),'Subject');
            assert(height(T)==73&&isequal(string(T.Subject),subjects));
            V(:,b,c)=double(T.Value);
        end
    end
    low=V(:,1,1)-V(:,1,2);high=V(:,2,1)-V(:,2,2);
    parts={low,high,high-low};
    for q=1:3
        endpoint=endpoint+1;D(:,endpoint)=parts{q};Metric(endpoint)=metrics(m);
        Contrast(endpoint)=labels(q);Unit(endpoint)=units(m);
    end
end
FeatureStats=runFamily(D,Contrast,Unit,2,p.Results.nPermutations, ...
    p.Results.nBootstrap,p.Results.minimumPairedN,p.Results.seed+100);
FeatureStats=addvars(FeatureStats,Metric,'Before',1,'NewVariableNames','Metric');

outputRoot=char(string(p.Results.outputRoot));
if isempty(outputRoot),outputRoot=char(cfgBurst.runDir);end
timeTag=char(datetime('now','Format','yyyyMMdd_HHmmss'));
runID=['BURST_BETA_SUBBAND_DURATION_MATCHED_' timeTag];
runDir=fullfile(outputRoot,runID);assert(~isfolder(runDir),'Refusing overwrite.');
[ok,msg]=mkdir(runDir);assert(ok,'Cannot create output: %s',msg);
cfgMatched=struct('runID',runID,'runDir',string(runDir), ...
    'burstFullFile',string(burstFullFile),'matchedMinimumSamples',matchedSamples, ...
    'matchedMinimumSeconds',matchedSamples/double(cfgBurst.samplingRateHz), ...
    'primaryThreshold',2,'sensitivityThresholds',[1.5 2.5], ...
    'nPermutations',p.Results.nPermutations, ...
    'nBootstrap',p.Results.nBootstrap,'seed',p.Results.seed, ...
    'analysisStatus','METHOD_MOTIVATED_ABSOLUTE_DURATION_SENSITIVITY');
DurationMatched=struct('MatchQC',MatchQC,'PrevalenceStats',PrevalenceStats, ...
    'FeatureStats',FeatureStats,'SubjectFeatures',SubjectFeatures);
save(fullfile(runDir,[runID '_RESULTS.mat']),'DurationMatched','cfgMatched','-v7.3');
writetable(MatchQC,fullfile(runDir,'DURATION_MATCH_QC.csv'));
writetable(PrevalenceStats,fullfile(runDir,'PREVALENCE_STATS.csv'));
writetable(FeatureStats,fullfile(runDir,'FEATURE_STATS.csv'));
writetable(SubjectFeatures,fullfile(runDir,'SUBJECT_FEATURES.csv'));
fprintf('\nBeta-subband duration-matched sensitivity complete:\n%s\n',runDir);
fprintf('\nMATCH QC\n');disp(MatchQC)
fprintf('\nPREVALENCE\n');disp(PrevalenceStats)
fprintf('\nFEATURES\n');disp(FeatureStats)
end

function keys=makeKeys(T)
keys=string(T.Subject)+"|"+string(T.Signal)+"|"+string(T.Band)+"|"+ ...
    string(T.Condition)+"|"+string(T.OriginalTrial)+"|"+string(T.ZThreshold);
end

function T=runFamily(D,contrast,unit,zThreshold,nPerm,nBoot,minN,seed)
k=size(D,2);N=sum(isfinite(D),1)';eligible=N>=minN;
mu=mean(D,1,'omitnan');sd=std(D,0,1,'omitnan');tObs=zeros(1,k);
for j=find(eligible)'
    if sd(j)>0,tObs(j)=mu(j)/(sd(j)/sqrt(N(j)));end
end
rng(double(seed),'twister');raw=zeros(1,k);adj=raw;
for first=1:1000:nPerm
    q=min(1000,nPerm-first+1);signs=2*(rand(q,size(D,1))>.5)-1;
    permT=zeros(q,k);
    for j=find(eligible)'
        ok=isfinite(D(:,j));X=signs(:,ok).*D(ok,j)';
        nm=mean(X,2);ns=std(X,0,2);permT(:,j)=nm./(ns/sqrt(sum(ok)));
        permT(~isfinite(permT(:,j)),j)=0;
        raw(j)=raw(j)+sum(abs(permT(:,j))>=abs(tObs(j)));
    end
    maxAbsT=max(abs(permT(:,eligible)),[],2);
    for j=find(eligible)',adj(j)=adj(j)+sum(maxAbsT>=abs(tObs(j)));end
end
rawP=nan(k,1);maxP=rawP;rawP(eligible)=(raw(eligible)+1)/(nPerm+1);
maxP(eligible)=(adj(eligible)+1)/(nPerm+1);dz=(mu./sd)';dz(~isfinite(dz))=NaN;
rng(double(seed)+1000,'twister');CI=nan(k,2);
for j=1:k
    d=D(isfinite(D(:,j)),j);n=numel(d);if n==0,continue,end
    boot=zeros(nBoot,1);done=0;
    while done<nBoot
        q=min(2000,nBoot-done);idx=randi(n,n,q);
        boot(done+(1:q))=mean(reshape(d(idx),n,q),1)';done=done+q;
    end
    CI(j,:)=prctile(boot,[2.5 97.5]);
end
proportionPositive=(sum(D>0,1)./N')';
T=table(contrast(:),unit(:),repmat(zThreshold,k,1),N,eligible,mu', ...
    CI(:,1),CI(:,2),dz,rawP,maxP,proportionPositive, ...
    'VariableNames',{'Contrast','Unit','ZThreshold','NPaired', ...
    'InferentiallyEligible','MeanContrast','BootstrapCI95Lower', ...
    'BootstrapCI95Upper','CohenDz','RawP','JointMaxTFWERP', ...
    'ProportionSubjectsPositive'});
T.SignificantJointMaxT05=T.InferentiallyEligible&T.JointMaxTFWERP<.05;
end
