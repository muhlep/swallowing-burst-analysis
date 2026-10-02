function [FeatureStats,SubjectFeatures,cfgFeatures] = burst_event_features(BurstFull,varargin)
% Subject-level cortical burst features for theta, alpha and beta at z=2.
% Count uses all valid ROI-trials. Duration, PeakZ and AUC use only
% non-left- and non-right-censored events. Aggregation: within ROI first,
% then median across ROIs. Max-|T| FWER is joint across eligible endpoints.

p=inputParser;
p.addParameter('outputRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('runLabel','EVENT_FEATURES',@(x)ischar(x)||isstring(x));
p.addParameter('nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
p.addParameter('nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
p.addParameter('minimumPairedN',20,@(x)isnumeric(x)&&isscalar(x)&&x>=10);
p.addParameter('randomSeed',20260817,@(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});

I=BurstFull.CorticalInventory; E=BurstFull.CorticalEvents;
subjects=sort(unique(string(I.Subject))); bands=["theta";"alpha";"beta"];
conditions=["EMG2";"EMG3"]; metrics=["BurstCountPerROITrial";"DurationSec";"PeakZ";"AUCAboveThreshold"];
assert(numel(subjects)==73,'Expected 73 subjects.');
assert(all(ismember({'Subject','Signal','Band','Condition','ZThreshold','TrialValid','NBurstsPrimary'},I.Properties.VariableNames)),'Inventory fields missing.');
assert(all(ismember({'Subject','Signal','Band','Condition','ZThreshold','DurationSec','PeakZ','AUCAboveThreshold','LeftCensored','RightCensored'},E.Properties.VariableNames)),'Event fields missing.');

rows=cell(73*3*2*4,10); row=0;
for s=1:73
 sid=subjects(s);
 for b=1:3
  band=bands(b);
  for c=1:2
   cond=conditions(c);
   iq=string(I.Subject)==sid & string(I.Band)==band & string(I.Condition)==cond & I.ZThreshold==2 & I.TrialValid;
   inv=I(iq,:); assert(height(inv)>0,'No inventory: %s %s %s',sid,band,cond);
   roiNames=unique(string(inv.Signal)); roiCount=nan(numel(roiNames),1);
   for r=1:numel(roiNames)
    roiCount(r)=mean(double(inv.NBurstsPrimary(string(inv.Signal)==roiNames(r))));
   end
   row=row+1; rows(row,:)={sid,cond,band,metrics(1),median(roiCount,'omitnan'),numel(roiCount),height(inv),sum(inv.NBurstsPrimary),false,"ALL_VALID_ROI_TRIALS"};

   eq=string(E.Subject)==sid & string(E.Band)==band & string(E.Condition)==cond & E.ZThreshold==2 & ~E.LeftCensored & ~E.RightCensored;
   ev=E(eq,:);
   fields={'DurationSec','PeakZ','AUCAboveThreshold'};
   for m=1:3
    roiValue=nan(numel(roiNames),1); nROIWithEvents=0;
    for r=1:numel(roiNames)
     q=string(ev.Signal)==roiNames(r);
     if any(q)
      nROIWithEvents=nROIWithEvents+1;
      roiValue(r)=median(double(ev.(fields{m})(q)),'omitnan');
     end
    end
    value=median(roiValue,'omitnan');
    row=row+1; rows(row,:)={sid,cond,band,metrics(m+1),value,nROIWithEvents,height(inv),height(ev),true,"COMPLETE_EVENTS_ROI_FIRST"};
   end
  end
 end
end
SubjectFeatures=cell2table(rows,'VariableNames',{'Subject','Condition','Band','Metric','Value', ...
 'NROIsContributing','NValidROITrials','NObservations','ConditionalOnBurst','AggregationRule'});

nEndpoint=12; D=nan(73,nEndpoint); V2=D; V3=D; Band=strings(nEndpoint,1); Metric=Band; endpoint=0;
for b=1:3
 for m=1:4
  endpoint=endpoint+1; Band(endpoint)=bands(b); Metric(endpoint)=metrics(m);
  A=sortrows(SubjectFeatures(string(SubjectFeatures.Band)==bands(b) & string(SubjectFeatures.Metric)==metrics(m) & string(SubjectFeatures.Condition)=="EMG2",:),'Subject');
  B=sortrows(SubjectFeatures(string(SubjectFeatures.Band)==bands(b) & string(SubjectFeatures.Metric)==metrics(m) & string(SubjectFeatures.Condition)=="EMG3",:),'Subject');
  assert(height(A)==73&&height(B)==73&&isequal(string(A.Subject),subjects)&&isequal(string(B.Subject),subjects),'Pairing failed.');
  V2(:,endpoint)=double(A.Value); V3(:,endpoint)=double(B.Value); D(:,endpoint)=V2(:,endpoint)-V3(:,endpoint);
 end
end

rng(double(p.Results.randomSeed),'twister');
[R,CI,Npaired,eligible]=maxTMissing(D,double(p.Results.nPermutations),double(p.Results.nBootstrap),double(p.Results.minimumPairedN));
FeatureStats=table(Band,Metric,Npaired,eligible,mean(V2,1,'omitnan')',mean(V3,1,'omitnan')', ...
 mean(D,1,'omitnan')',CI(:,1),CI(:,2),R.Dz,R.RawP,R.MaxTP, ...
 'VariableNames',{'Band','Metric','NPaired','InferentiallyEligible','MeanEMG2', ...
 'MeanEMG3','MeanPairedDifference','BootstrapCI95Lower','BootstrapCI95Upper', ...
 'PairedCohensDz','PermutationP','JointMaxTFWERP'});
FeatureStats.SignificantJointMaxT05=FeatureStats.InferentiallyEligible & FeatureStats.JointMaxTFWERP<.05;
FeatureStats.Interpretation=repmat("INFERENTIAL",nEndpoint,1);
FeatureStats.Interpretation(~eligible)="DESCRIPTIVE_TOO_FEW_PAIRED_SUBJECTS";

outputRoot=char(string(p.Results.outputRoot)); if isempty(outputRoot), outputRoot=pwd; end
runID=sprintf('BURST_%s_%s',char(string(p.Results.runLabel)),char(datetime('now','Format','yyyyMMdd_HHmmss')));
runDir=fullfile(outputRoot,runID); assert(~isfolder(runDir),'Refusing overwrite: %s',runDir);
[ok,msg]=mkdir(runDir); assert(ok,'Cannot create output: %s',msg);
cfgFeatures=struct('runID',runID,'runDir',runDir,'bands',bands,'threshold',2, ...
 'metrics',metrics,'minimumPairedN',p.Results.minimumPairedN, ...
 'nPermutations',p.Results.nPermutations,'nBootstrap',p.Results.nBootstrap, ...
 'randomSeed',p.Results.randomSeed,'multiplicity', ...
 'joint max-|T| FWER across all inferentially eligible standard-band feature endpoints', ...
 'censoring','count includes right-censored primary bursts; conditional features exclude left and right censoring', ...
 'aggregation','subject x condition x band: ROI summary first, then median across ROIs');
save(fullfile(runDir,[runID '_STANDARD_BURST_FEATURES.mat']), ...
 'FeatureStats','SubjectFeatures','cfgFeatures','-v7.3');
writetable(FeatureStats,fullfile(runDir,[runID '_FEATURE_STATS.csv']));
writetable(SubjectFeatures,fullfile(runDir,[runID '_SUBJECT_FEATURES.csv']));
fprintf('\nEvent feature analysis complete: %s\n',runDir); disp(FeatureStats);
end

function [O,CI,N,eligible]=maxTMissing(D,nPerm,nBoot,minN)
k=size(D,2); N=sum(isfinite(D),1)'; eligible=N>=minN; tObs=zeros(1,k); mu=nan(1,k); sd=nan(1,k);
for j=1:k
 d=D(isfinite(D(:,j)),j); mu(j)=mean(d); sd(j)=std(d); if eligible(j)&&sd(j)>0, tObs(j)=mu(j)/(sd(j)/sqrt(numel(d))); end
end
raw=zeros(1,k); adj=zeros(1,k);
for first=1:1000:nPerm
 q=min(1000,nPerm-first+1); signs=2*(rand(q,size(D,1))>.5)-1; TP=zeros(q,k);
 for j=find(eligible)'
  ok=isfinite(D(:,j)); X=signs(:,ok).*D(ok,j)'; nm=mean(X,2); ns=std(X,0,2); TP(:,j)=nm./(ns/sqrt(sum(ok))); TP(~isfinite(TP(:,j)),j)=0;
  raw(j)=raw(j)+sum(abs(TP(:,j))>=abs(tObs(j)));
 end
 mx=max(abs(TP(:,eligible)),[],2);
 for j=find(eligible)', adj(j)=adj(j)+sum(mx>=abs(tObs(j))); end
end
rawP=nan(k,1); maxP=rawP; rawP(eligible)=(raw(eligible)+1)/(nPerm+1); maxP(eligible)=(adj(eligible)+1)/(nPerm+1);
dz=(mu./sd)'; dz(~isfinite(dz))=NaN; O=struct('Dz',dz,'RawP',rawP,'MaxTP',maxP); CI=nan(k,2);
for j=1:k
 d=D(isfinite(D(:,j)),j); n=numel(d); if n==0, continue; end
 boot=zeros(nBoot,1); done=0;
 while done<nBoot
  q=min(2000,nBoot-done); idx=randi(n,n,q); sampled=reshape(d(idx),n,q);
  boot(done+(1:q))=mean(sampled,1)'; done=done+q;
 end
 CI(j,:)=prctile(boot,[2.5 97.5]);
end
end
