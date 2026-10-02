function [ROIPilot,cfgROI] = burst_extract_roi_timeseries(PilotCommonFilter,cfgCommonFilterPilot,varargin)
% Extract one sign-deterministic ROI PC per ROI using weights learned jointly
% from retained EMG2 and EMG3 trials. The same weights are applied to both.

p=inputParser;
p.addParameter('runLabel','ROI_TIMESERIES',@(x)ischar(x)||isstring(x));
p.parse(varargin{:});
assert(isstruct(PilotCommonFilter)&&isstruct(cfgCommonFilterPilot),'Invalid pilot inputs.');
assert(isfield(PilotCommonFilter,'sourceCondition')&&numel(PilotCommonFilter.sourceCondition)==2, ...
 'Two source conditions are required.');
assert(isfield(PilotCommonFilter,'ROIs')&&numel(PilotCommonFilter.ROIs)==21,'Expected 21 ROIs.');

conditions=["EMG2";"EMG3"]; sources=PilotCommonFilter.sourceCondition;
ROIs=PilotCommonFilter.ROIs; nROI=numel(ROIs); nTime=600;
runID=sprintf('BURST_%s_%s',char(string(p.Results.runLabel)), ...
 char(datetime('now','Format','yyyyMMdd_HHmmss')));
runDir=fullfile(cfgCommonFilterPilot.runDir,runID);
assert(~isfolder(runDir),'Refusing overwrite: %s',runDir);
[ok,msg]=mkdir(runDir); assert(ok,'Cannot create %s: %s',runDir,msg);

nTrials=[numel(sources{1}.trial) numel(sources{2}.trial)];
roiData={nan(nROI,nTrials(1),nTime),nan(nROI,nTrials(2),nTime)};
roiWeights=cell(nROI,1); roiCenters=cell(nROI,1); roiIndices=cell(nROI,1);
QC=cell(nROI,10);

for r=1:nROI
 idx=double(ROIs(r).idx(:)); idx=idx(isfinite(idx)&idx>=1);
 assert(~isempty(idx),'ROI %s has no voxel indices.',string(ROIs(r).name));
 Xparts=cell(sum(nTrials),1); q=0;
 for ci=1:2
  for t=1:nTrials(ci)
   q=q+1; X=nan(numel(idx),nTime);
   for v=1:numel(idx)
    assert(idx(v)<=numel(sources{ci}.trial(t).mom),'ROI index outside moment cell.');
    x=double(sources{ci}.trial(t).mom{idx(v)});
    assert(~isempty(x)&&numel(x)==nTime,'Missing/invalid moment: ROI %s, %s, trial %d.', ...
     string(ROIs(r).name),conditions(ci),t);
    X(v,:)=reshape(x,1,[]);
   end
   assert(all(isfinite(X),'all'),'Non-finite ROI moments.'); Xparts{q}=X;
  end
 end
 Xall=horzcat(Xparts{:}); center=mean(Xall,2); Xcenter=Xall-center;
 [U,S,~]=svd(Xcenter,'econ'); w=U(:,1);
 [~,anchor]=max(abs(w)); if w(anchor)<0, w=-w; end
 latent=diag(S).^2; explained=100*latent(1)/sum(latent);

 for ci=1:2
  for t=1:nTrials(ci)
   X=nan(numel(idx),nTime);
   for v=1:numel(idx), X(v,:)=reshape(double(sources{ci}.trial(t).mom{idx(v)}),1,[]); end
   roiData{ci}(r,t,:)=reshape(w'*(X-center),1,1,nTime);
  end
 end

 C=corr(Xcenter','Rows','pairwise');
 if numel(idx)>1
  upper=C(triu(true(size(C)),1)); medAbsCorr=median(abs(upper),'omitnan');
 else
  medAbsCorr=1;
 end
 rms2=sqrt(mean(reshape(roiData{1}(r,:,:),[],1).^2,'omitnan'));
 rms3=sqrt(mean(reshape(roiData{2}(r,:,:),[],1).^2,'omitnan'));
 roiWeights{r}=w; roiCenters{r}=center; roiIndices{r}=idx;
 QC(r,:)={string(ROIs(r).name),numel(idx),explained,medAbsCorr,anchor, ...
  min(w),max(w),rms2,rms3,rms2/rms3};
 fprintf('%02d/%02d | %-12s | vox=%d | PC1=%.2f%%\n',r,nROI,string(ROIs(r).name),numel(idx),explained);
end

ROIPCAQC=cell2table(QC,'VariableNames',{'ROI','NVoxels','PC1ExplainedPercent', ...
 'MedianAbsoluteIntervoxelCorrelation','AnchorVoxelWithinROI','MinimumWeight', ...
 'MaximumWeight','RMSEMG2','RMSEMG3','RMSEMG2toEMG3Ratio'});
assert(all(isfinite(roiData{1}),'all')&&all(isfinite(roiData{2}),'all'),'Non-finite ROI output.');
assert(isequal(size(roiData{1}),[21 nTrials(1) 600]),'Unexpected EMG2 ROI dimensions.');
assert(isequal(size(roiData{2}),[21 nTrials(2) 600]),'Unexpected EMG3 ROI dimensions.');

trialInventory=PilotCommonFilter.TrialInventory;
time={double(sources{1}.time(:)'),double(sources{2}.time(:)')};
roiLabels=string({ROIs.name})';
emgData=cell(2,1);
emgEnvelope=cell(2,1);
for ci=1:2
 assert(isfield(PilotCommonFilter,'dataEMG')&&numel(PilotCommonFilter.dataEMG)>=ci, ...
  'Compact EMG data missing for condition %s.',conditions(ci));
 emgData{ci}=nan(1,nTrials(ci),nTime);
 emgEnvelope{ci}=nan(1,nTrials(ci),nTime);
 for t=1:nTrials(ci)
  x=double(PilotCommonFilter.dataEMG{ci}.trial{t});
  assert(size(x,1)==1&&size(x,2)==nTime,'Unexpected EMG dimensions.');
  emgData{ci}(1,t,:)=reshape(x,1,1,nTime);
  assert(isfield(PilotCommonFilter,'dataEMGEnvelope'), ...
   'Padded EMG envelope is missing. Update burst_source_reconstruction.');
  e=double(PilotCommonFilter.dataEMGEnvelope{ci}.trial{t});
  assert(size(e,1)==1&&size(e,2)==nTime&&all(isfinite(e),'all'), ...
   'Unexpected EMG-envelope dimensions or values.');
  emgEnvelope{ci}(1,t,:)=reshape(e,1,1,nTime);
 end
end
ROIPilot=struct('subject',cfgCommonFilterPilot.subject,'band',cfgCommonFilterPilot.band, ...
 'bandLabel',cfgCommonFilterPilot.bandLabel,'conditions',conditions,'roiData',{roiData}, ...
 'emgRawTarget',{emgData},'emgProcessing','raw unfiltered target-window signal', ...
 'emgEnvelopeTarget',{emgEnvelope}, ...
 'emgEnvelopeProcessing','20-200 Hz, full-wave rectified, 10-Hz low-pass; computed with 1-s padding', ...
 'roiLabels',roiLabels,'roiIndices',{roiIndices},'roiWeights',{roiWeights}, ...
 'roiCenters',{roiCenters},'time',{time},'trialInventory',trialInventory, ...
 'ROIPCAQC',ROIPCAQC,'weightDefinition', ...
 'joint EMG2+EMG3 PC1; row-centered; anchor voxel weight forced positive');
cfgROI=struct('runID',runID,'runDir',runDir,'parentRunID',cfgCommonFilterPilot.runID, ...
 'method','joint_sign_deterministic_roi_pc1','nTime',nTime,'nTrials',nTrials);
save(fullfile(runDir,[runID '_ROI_PCA_PILOT.mat']),'ROIPilot','cfgROI','-v7.3');
writetable(ROIPCAQC,fullfile(runDir,[runID '_ROI_PCA_QC.csv']));
writetable(trialInventory,fullfile(runDir,[runID '_TRIAL_INVENTORY.csv']));
fprintf('\nROI time-series extraction complete: %s\n',runDir); disp(ROIPCAQC);
end
