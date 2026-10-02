function [EMGROIStats,cfgEMG] = burst_cortical_emg_overlap(BurstFull,varargin)
% ROI-specific cortical overlap with the swallowing EMG envelope at z=2.
% Run once with theta/alpha/beta and once with low/high gamma so that the
% two prespecified correction families remain separate.
% Exact within-subject null means:
%   trial shuffle = average over all cortical-trial x EMG-trial pairings
%   circular shift = average over all non-zero within-trial EMG shifts
% Group inference uses subject sign-flip max-|T| across every band,
% condition, and ROI in the selected family, separately for both nulls.

p=inputParser;
p.addParameter('outputRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('bands',["theta";"alpha";"beta"], ...
    @(x)isstring(x)||iscellstr(x)||ischar(x));
p.addParameter('runLabel','CORTICAL_EMG_OVERLAP',@(x)ischar(x)||isstring(x));
p.addParameter('nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
p.addParameter('nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
p.addParameter('randomSeed',20260817,@(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
CI=BurstFull.CorticalInventory; CE=BurstFull.CorticalEvents;
MI=BurstFull.EMGInventory; ME=BurstFull.EMGEvents;
subjects=sort(unique(string(CI.Subject))); rois=sort(unique(string(CI.Signal)));
bands=string(p.Results.bands(:)); conditions=["EMG2";"EMG3"];
validStandard=isequal(bands,["theta";"alpha";"beta"]);
validGamma=isequal(bands,["lowgamma";"highgamma"]);
assert(validStandard||validGamma, ...
    'Bands must define either theta/alpha/beta or low/high gamma.');
assert(numel(subjects)==73&&numel(rois)==21,'Expected 73 subjects and 21 ROIs.');
threshold=2; nSamples=600; nEndpoint=numel(bands)*2*21;
R=cell(73*nEndpoint,17); row=0; timer=tic;

for s=1:73
 sid=subjects(s);
 for b=1:numel(bands)
  band=bands(b);
  for c=1:2
   cond=conditions(c);
   iq=string(CI.Subject)==sid&string(CI.Band)==band&string(CI.Condition)==cond&CI.ZThreshold==threshold;
   mq=string(MI.Subject)==sid&string(MI.Condition)==cond&MI.ZThreshold==threshold;
   trials=intersect(unique(double(CI.OriginalTrial(iq))),unique(double(MI.OriginalTrial(mq))));
   nTrial=numel(trials); assert(nTrial>0,'No common trials: %s | %s.',sid,cond);
   Mmask=false(nTrial,nSamples);
   emq=string(ME.Subject)==sid&string(ME.Condition)==cond&ME.ZThreshold==threshold&ME.PrimaryEligible;
   em=ME(emq,:); [keepM,tiM]=ismember(double(em.OriginalTrial),trials); em=em(keepM,:); tiM=tiM(keepM);
   for k=1:height(em)
    a=max(1,double(em.OnsetSample(k))); z=min(nSamples,double(em.OffsetSample(k)));
    Mmask(tiM(k),a:z)=true;
   end
   emgOcc=mean(Mmask,'all');
   egq=string(CE.Subject)==sid&string(CE.Band)==band&string(CE.Condition)==cond& ...
    CE.ZThreshold==threshold&CE.PrimaryEligible;
   corticalGroupEvents=CE(egq,:);

   for r=1:21
    roi=rois(r); Cmask=false(nTrial,nSamples);
    ev=corticalGroupEvents(string(corticalGroupEvents.Signal)==roi,:);
    [keepC,tiC]=ismember(double(ev.OriginalTrial),trials); ev=ev(keepC,:); tiC=tiC(keepC);
    for k=1:height(ev)
     a=max(1,double(ev.OnsetSample(k))); z=min(nSamples,double(ev.OffsetSample(k)));
     Cmask(tiC(k),a:z)=true;
    end
    corticalOcc=mean(Cmask,'all'); observed=mean(Cmask&Mmask,'all');

    pairOverlap=(double(Cmask)*double(Mmask)')/nSamples;
    trialNull=mean(pairOverlap,'all');
    byShift=zeros(nTrial,nSamples);
    for t=1:nTrial
     byShift(t,:)=real(ifft(conj(fft(double(Cmask(t,:)))).*fft(double(Mmask(t,:)))))/nSamples;
    end
    assert(abs(mean(byShift(:,1))-observed)<1e-10,'Circular overlap QC failed.');
    circularNull=mean(byShift(:,2:end),'all');
    row=row+1;
    R(row,:)={sid,band,cond,roi,nTrial,threshold,100*corticalOcc,100*emgOcc, ...
     100*observed,100*trialNull,100*(observed-trialNull),100*circularNull, ...
     100*(observed-circularNull),sum(any(Cmask,2)),sum(any(Mmask,2)), ...
     100*corticalOcc*emgOcc,21};
   end
  end
 end
 fprintf('%d/73 | %s | elapsed %s\n',s,sid,char(duration(0,0,toc(timer),'Format','hh:mm:ss')));
end

SubjectROIOverlap=cell2table(R,'VariableNames',{'Subject','Band','Condition','ROI', ...
 'NTrials','ZThreshold','CorticalOccupancyPercent','EMGOccupancyPercent', ...
 'ObservedOverlapPercent','ExactTrialShuffleMeanPercent','TrialSpecificExcessPercent', ...
 'ExactCircularShiftMeanPercent','TimeLockedExcessPercent','NTrialsWithCorticalBurst', ...
 'NTrialsWithEMGBurst','IndependenceExpectedOverlapPercent','NetworkNROI'});

% Form participant-by-endpoint matrices in a fixed order.
DT=nan(73,nEndpoint); DC=DT; Obs=DT; CortOcc=DT; EMGOcc=DT;
EndpointBand=strings(nEndpoint,1); EndpointCondition=strings(nEndpoint,1); EndpointROI=strings(nEndpoint,1); e=0;
for b=1:numel(bands)
 for c=1:2
  for r=1:21
   e=e+1; q=string(SubjectROIOverlap.Band)==bands(b)& ...
    string(SubjectROIOverlap.Condition)==conditions(c)&string(SubjectROIOverlap.ROI)==rois(r);
   X=sortrows(SubjectROIOverlap(q,:),'Subject'); assert(height(X)==73&&isequal(string(X.Subject),subjects),'Endpoint pairing failed.');
   DT(:,e)=double(X.TrialSpecificExcessPercent); DC(:,e)=double(X.TimeLockedExcessPercent);
   Obs(:,e)=double(X.ObservedOverlapPercent); CortOcc(:,e)=double(X.CorticalOccupancyPercent);
   EMGOcc(:,e)=double(X.EMGOccupancyPercent); EndpointBand(e)=bands(b); EndpointCondition(e)=conditions(c); EndpointROI(e)=rois(r);
  end
 end
end
rng(double(p.Results.randomSeed),'twister');
[trialResult,trialCI]=maxTInference(DT,double(p.Results.nPermutations),double(p.Results.nBootstrap));
[timeResult,timeCI]=maxTInference(DC,double(p.Results.nPermutations),double(p.Results.nBootstrap));
EMGROIStats=table(EndpointBand,EndpointCondition,EndpointROI,repmat(73,nEndpoint,1), ...
 mean(CortOcc,1)',mean(EMGOcc,1)',mean(Obs,1)',mean(DT,1)',trialCI(:,1),trialCI(:,2), ...
 trialResult.Dz,trialResult.RawP,trialResult.MaxTP,mean(DC,1)',timeCI(:,1),timeCI(:,2), ...
 timeResult.Dz,timeResult.RawP,timeResult.MaxTP, ...
 'VariableNames',{'Band','Condition','ROI','N','MeanCorticalOccupancyPercent', ...
 'MeanEMGOccupancyPercent','MeanObservedOverlapPercent','MeanTrialSpecificExcessPercent', ...
 'TrialExcessCI95Lower','TrialExcessCI95Upper','TrialExcessDz','TrialExcessPermutationP', ...
 'TrialExcessMaxTFWERP','MeanTimeLockedExcessPercent','TimeLockedCI95Lower', ...
 'TimeLockedCI95Upper','TimeLockedDz','TimeLockedPermutationP','TimeLockedMaxTFWERP'});
EMGROIStats.TrialSpecificMaxT05=EMGROIStats.TrialExcessMaxTFWERP<.05;
EMGROIStats.TimeLockedMaxT05=EMGROIStats.TimeLockedMaxTFWERP<.05;
EMGROIStats.BothNullModelsMaxT05=EMGROIStats.TrialSpecificMaxT05&EMGROIStats.TimeLockedMaxT05;
Summary=groupsummary(EMGROIStats,{'Band','Condition'},'sum', ...
 {'TrialSpecificMaxT05','TimeLockedMaxT05','BothNullModelsMaxT05'});

outputRoot=char(string(p.Results.outputRoot)); if isempty(outputRoot), outputRoot=pwd; end
runID=sprintf('BURST_%s_%s',char(string(p.Results.runLabel)),char(datetime('now','Format','yyyyMMdd_HHmmss')));
runDir=fullfile(outputRoot,runID); assert(~isfolder(runDir),'Refusing overwrite: %s',runDir);
[ok,msg]=mkdir(runDir); assert(ok,'Cannot create output: %s',msg);
cfgEMG=struct('runID',runID,'runDir',runDir,'threshold',2,'nSubjects',73, ...
 'bands',bands,'nEndpoints',nEndpoint,'nPermutations',p.Results.nPermutations, ...
 'nBootstrap',p.Results.nBootstrap,'randomSeed',p.Results.randomSeed, ...
 'multiplicity',sprintf('separate max-|T| FWER across all %d endpoints for each null model',nEndpoint), ...
 'trialNull','exact mean over all within-subject/condition trial pairings', ...
 'timeNull','exact mean over all non-zero circular shifts');
save(fullfile(runDir,[runID '_RESULTS.mat']), ...
 'EMGROIStats','SubjectROIOverlap','Summary','cfgEMG','-v7.3');
writetable(EMGROIStats,fullfile(runDir,[runID '_ROI_STATS.csv']));
writetable(SubjectROIOverlap,fullfile(runDir,[runID '_SUBJECT_ROI_OVERLAP.csv']));
writetable(Summary,fullfile(runDir,[runID '_SUMMARY.csv']));
fprintf('\nCortical-EMG overlap analysis complete: %s\n',runDir); disp(Summary);
end

function [O,CI]=maxTInference(D,nPerm,nBoot)
n=size(D,1); mEndpoint=size(D,2); mu=mean(D,1); sd=std(D,0,1); tObs=zeros(1,mEndpoint);
valid=sd>0&isfinite(sd); tObs(valid)=mu(valid)./(sd(valid)/sqrt(n)); sumSq=sum(D.^2,1);
raw=zeros(1,mEndpoint); adj=raw; chunk=1000;
for first=1:chunk:nPerm
 m=min(chunk,nPerm-first+1); signs=2*(rand(m,n)>.5)-1; nm=(signs*D)/n;
 nv=(repmat(sumSq,m,1)-n*nm.^2)/(n-1); nt=zeros(size(nm)); ok=nv>0;
 nt(ok)=nm(ok)./sqrt(nv(ok)/n); maxAbs=max(abs(nt),[],2);
 raw=raw+sum(abs(nt)>=abs(tObs),1); adj=adj+sum(maxAbs>=abs(tObs),1);
end
O=struct('Dz',(mu./sd)','RawP',((raw+1)/(nPerm+1))','MaxTP',((adj+1)/(nPerm+1))');
O.Dz(~isfinite(O.Dz))=0; CI=nan(mEndpoint,2);
for j=1:mEndpoint
 d=D(:,j); boot=zeros(nBoot,1); done=0;
 while done<nBoot
  m=min(2000,nBoot-done); idx=randi(n,n,m); sampled=reshape(d(idx),n,m);
  boot(done+(1:m))=mean(sampled,1)'; done=done+m;
 end
 CI(j,:)=prctile(boot,[2.5 97.5]);
end
end
