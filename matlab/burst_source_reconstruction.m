function [Pilot,cfgPilot] = burst_source_reconstruction(QCLocked,LockedSubjectQC,EMGChannelMapping,TimeWindowDefinition,varargin)
% Single-participant/single-band common-filter reconstruction.
% Study-specific data locations are supplied by the caller.

p=inputParser;
p.addParameter('subject','',@(x)ischar(x)||isstring(x));
p.addParameter('band',[13 30],@(x)isnumeric(x)&&numel(x)==2);
p.addParameter('bandLabel','beta',@(x)ischar(x)||isstring(x));
p.addParameter('lambda','5%',@(x)ischar(x)||isstring(x)||isnumeric(x));
p.addParameter('originalCodeRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('outputRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('runLabel','SOURCE_RECONSTRUCTION',@(x)ischar(x)||isstring(x));
p.addParameter('saveFullPilot',true,@(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});

sid=string(p.Results.subject); band=double(p.Results.band(:)');
bandLabel=string(p.Results.bandLabel); conditions=["EMG2";"EMG3"];
assert(strlength(sid)>0,'A participant pairing key must be supplied.');
assert(isfolder(char(string(p.Results.originalCodeRoot))), ...
 'Study-specific configuration code directory is missing.');
assert(isfolder(char(string(p.Results.outputRoot))), ...
 'An existing output directory must be supplied.');
assert(istable(QCLocked)&&istable(LockedSubjectQC),'Locked QC inputs must be tables.');
assert(istable(EMGChannelMapping)&&istable(TimeWindowDefinition),'Locked mapping inputs must be tables.');
sr=string(LockedSubjectQC.Subject)==sid; assert(sum(sr)==1,'Subject not unique in locked QC.');
assert(~LockedSubjectQC.ExcludeSubject(sr),'Subject %s is excluded.',sid);

codeRoot=char(string(p.Results.originalCodeRoot)); addpath(codeRoot);
assert(exist('pmd15_cfg','file')==2,'pmd15_cfg not found.');
ftDefaultsFile=which('ft_defaults');
assert(~isempty(ftDefaultsFile),'ft_defaults is not on the MATLAB path.');
ftRoot=fileparts(ftDefaultsFile);
ctfDir=fullfile(ftRoot,'external','ctf');
assert(isfolder(ctfDir),'FieldTrip CTF reader directory missing: %s',ctfDir);
addpath(ctfDir); rehash toolboxcache;
assert(exist('readCTFds','file')==2,'readCTFds is unavailable after adding %s.',ctfDir);
c=pmd15_cfg('subj',char(sid));
mr=string(EMGChannelMapping.Subject)==sid; assert(sum(mr)==1,'EMG mapping not unique.');
dataset=char(EMGChannelMapping.Dataset(mr)); emgLabel=char(EMGChannelMapping.EMGChannel(mr));
assert(isfolder(dataset),'Dataset missing: %s',dataset);

stamp=char(datetime('now','Format','yyyyMMdd_HHmmss'));
runID=sprintf('BURST_%s_%s',char(string(p.Results.runLabel)),stamp);
runDir=fullfile(char(string(p.Results.outputRoot)),runID);
assert(~isfolder(runDir),'Refusing overwrite: %s',runDir);
[ok,msg]=mkdir(runDir); assert(ok,'Cannot create %s: %s',runDir,msg);
cfgPilot=struct('runID',runID,'runDir',runDir,'subject',sid,'band',band, ...
 'bandLabel',bandLabel,'lambda',p.Results.lambda,'conditions',conditions, ...
 'filterDefinition','equal-weight mean of EMG2 and EMG3 condition covariance matrices', ...
 'unitStandard','mm','rejectionSource','locked sensor-level quality control');
save(fullfile(runDir,[runID '_LOCKED_INPUTS.mat']),'cfgPilot','QCLocked', ...
 'LockedSubjectQC','EMGChannelMapping','TimeWindowDefinition','-v7.3');

% Load geometry and normalize units in memory.
H=load(c.headmodel_file); vol=pickfield(H,{'vol','headmodel'});
G=load(c.warped_grid_file); grid=pickfield(G,{'warped_sourcemodel','sourcemodel'});
R=load(c.roi_file,'ROIs'); assert(isfield(R,'ROIs')&&numel(R.ROIs)==21,'Expected 21 ROIs.');
ROIs=R.ROIs; vol=ft_convert_units(vol,'mm'); grid=ft_convert_units(grid,'mm');
assert(strcmpi(vol.unit,'mm')&&strcmpi(grid.unit,'mm'),'Unit conversion failed.');

dataBand=cell(2,1); dataEMG=cell(2,1); TLavg=cell(2,1); TLtrial=cell(2,1);
dataBand=cell(2,1); dataEMG=cell(2,1); dataEMGEnvelope=cell(2,1); TLavg=cell(2,1); TLtrial=cell(2,1);
TrialInventory=cell(2,1); timer=tic;
for ci=1:2
 cond=conditions(ci); wr=string(TimeWindowDefinition.Condition)==cond;
 assert(sum(wr)==1,'Window definition missing for %s.',cond);
 target=[TimeWindowDefinition.WindowStartSec(wr) TimeWindowDefinition.WindowEndSec(wr)];
 loadWindow=target+[-1 1];
 fprintf('\n%s | %s | target [%g %g] | load [%g %g]\n',sid,cond,target,loadWindow);

 cfg=[]; cfg.dataset=dataset; cfg.continuous='yes'; cfg.trialfun='ft_trialfun_general';
 cfg.headerformat='ctf_ds'; cfg.dataformat='ctf_ds'; cfg.eventformat='ctf_ds';
 cfg.trialdef.eventtype=char(cond); cfg.trialdef.prestim=abs(loadWindow(1));
 cfg.trialdef.poststim=loadWindow(2); cfg=ft_definetrial(cfg); originalTrl=cfg.trl;
 cfg.channel={'MEG',emgLabel}; cfg.demean='no'; D=ft_preprocessing(cfg);
 nOriginal=numel(D.trial); D.trialinfo=(1:nOriginal)';

 qr=string(QCLocked.Subject)==sid & string(QCLocked.Condition)==cond;
 Qt=sortrows(QCLocked(qr,:), 'OriginalTrial');
 assert(height(Qt)==nOriginal,'QC/raw trial mismatch for %s: %d vs %d.',cond,height(Qt),nOriginal);
 assert(isequal(double(Qt.OriginalTrial),(1:nOriginal)'),'Original trial keys are not consecutive.');
 keep=~Qt.ExcludeTrialGeneralArtifact;
 cfgs=[]; cfgs.trials=find(keep); D=ft_selectdata(cfgs,D); D.trialinfo=find(keep);

 cfgf=[]; cfgf.channel='MEG'; cfgf.bpfilter='yes'; cfgf.bpfreq=band;
 cfgf.bpfilttype='but'; cfgf.bpfiltdir='twopass'; cfgf.demean='no';
 DB=ft_preprocessing(cfgf,D); DB=crophalfopen(DB,target);
 cfgemg=[]; cfgemg.channel={emgLabel}; DE=ft_selectdata(cfgemg,D); DE=crophalfopen(DE,target);
 cfgemg=[]; cfgemg.channel={emgLabel}; DEpad=ft_selectdata(cfgemg,D);
 DEenv=makeemgenvelope(DEpad); DEenv=crophalfopen(DEenv,target);
 DE=crophalfopen(DEpad,target);
 assert(all(cellfun(@(x)size(x,2)==600,DB.trial)),'MEG target is not 600 samples.');
 assert(all(cellfun(@(x)size(x,2)==600,DE.trial)),'EMG target is not 600 samples.');

 cfgt=[]; cfgt.covariance='yes'; cfgt.covariancewindow='all'; cfgt.keeptrials='no';
 TLavg{ci}=ft_timelockanalysis(cfgt,DB);
 cfgt.keeptrials='yes'; TLtrial{ci}=ft_timelockanalysis(cfgt,DB);
 dataBand{ci}=DB; dataEMG{ci}=DE;
 dataBand{ci}=DB; dataEMG{ci}=DE; dataEMGEnvelope{ci}=DEenv;
 TrialInventory{ci}=table(repmat(sid,sum(keep),1),repmat(cond,sum(keep),1), ...
  find(keep),originalTrl(keep,1)-originalTrl(keep,3), ...
  'VariableNames',{'Subject','Condition','OriginalTrial','TriggerSample'});
 fprintf('retained %d/%d | elapsed %s\n',sum(keep),nOriginal, ...
  char(duration(0,0,toc(timer),'Format','hh:mm:ss')));
end

assert(isequal(TLavg{1}.label,TLavg{2}.label),'MEG channel labels differ across conditions.');
commonTL=TLavg{1}; commonTL.cov=(double(TLavg{1}.cov)+double(TLavg{2}.cov))/2;
commonTL.avg=(double(TLavg{1}.avg)+double(TLavg{2}.avg))/2;
if isfield(TLavg{1},'var')&&isfield(TLavg{2},'var')
 commonTL.var=(double(TLavg{1}.var)+double(TLavg{2}.var))/2;
end

cfgsm=[]; cfgsm.sourcemodel=grid; cfgsm.headmodel=vol; cfgsm.grad=dataBand{1}.grad;
cfgsm.channel=cellstr(commonTL.label); sourcemodel=ft_prepare_sourcemodel(cfgsm);
cfglf=[]; cfglf.headmodel=vol; cfglf.sourcemodel=sourcemodel;
cfglf.grad=dataBand{1}.grad; cfglf.channel=cellstr(commonTL.label);
cfglf.reducerank=2; cfglf.normalize='yes';
leadfield=ft_prepare_leadfield(cfglf,commonTL);

cfgsa=[]; cfgsa.method='lcmv'; cfgsa.sourcemodel=leadfield; cfgsa.headmodel=vol;
cfgsa.lcmv.keepfilter='yes'; cfgsa.lcmv.fixedori='yes';
cfgsa.lcmv.projectmom='yes'; cfgsa.lcmv.lambda=p.Results.lambda;
sourceCommon=ft_sourceanalysis(cfgsa,commonTL);
filters=getfilters(sourceCommon); inside=find(sourceCommon.inside);
validInside=inside(~cellfun(@isempty,filters(inside)));
filterNorm=nan(numel(validInside),1);
for k=1:numel(validInside), filterNorm(k)=norm(double(filters{validInside(k)}(:))); end

% Apply the exact same precomputed filter to each condition.
sourceCondition=cell(2,1);
for ci=1:2
 applyLF=leadfield; applyLF.filter=filters;
 cfga=[]; cfga.method='lcmv'; cfga.sourcemodel=applyLF; cfga.headmodel=vol;
 cfga.keeptrials='yes'; cfga.rawtrial='yes'; cfga.lcmv.keepmom='yes';
 cfga.lcmv.projectmom='yes'; cfga.lcmv.fixedori='yes';
 sourceCondition{ci}=ft_sourceanalysis(cfga,TLtrial{ci});
end

roiVoxelCounts=arrayfun(@(x)numel(x.idx),ROIs(:));
roiValidFilterCounts=zeros(numel(ROIs),1);
for r=1:numel(ROIs)
 idx=ROIs(r).idx(:); idx=idx(idx>=1 & idx<=numel(filters));
 roiValidFilterCounts(r)=sum(~cellfun(@isempty,filters(idx)));
end
ROIQC=table(string({ROIs.name})',roiVoxelCounts,roiValidFilterCounts, ...
 'VariableNames',{'ROI','NVoxels','NValidFilters'});
FilterQC=table(sid,bandLabel,numel(inside),numel(validInside), ...
 min(filterNorm,[],'omitnan'),median(filterNorm,'omitnan'),max(filterNorm,[],'omitnan'), ...
 rcond(double(commonTL.cov)),numel(dataBand{1}.trial),numel(dataBand{2}.trial), ...
 'VariableNames',{'Subject','Band','NInside','NValidFilters','MinFilterNorm', ...
 'MedianFilterNorm','MaxFilterNorm','CovarianceRcond','NTrialsEMG2','NTrialsEMG3'});
assert(all(ROIQC.NValidFilters>0),'At least one ROI has no valid spatial filter.');
assert(FilterQC.NValidFilters==FilterQC.NInside,'Some inside-grid filters are missing.');
assert(all(isfinite(filterNorm)),'Non-finite filter norm.');

Pilot=struct('FilterQC',FilterQC,'ROIQC',ROIQC, ...
 'TrialInventory',vertcat(TrialInventory{:}), ...
 'commonTL',commonTL,'leadfield',leadfield, ...
 'sourceCommon',sourceCommon, ...
 'sourceCondition',{sourceCondition}, ...
 'dataEMG',{dataEMG}, ...
 'dataEMGEnvelope',{dataEMGEnvelope}, ...
 'ROIs',ROIs);

if p.Results.saveFullPilot
 save(fullfile(runDir, ...
  [runID '_COMMON_FILTER_PILOT.mat']), ...
  'Pilot','cfgPilot','-v7.3');
else
 TrialInventoryCompact=Pilot.TrialInventory;
 save(fullfile(runDir, ...
  [runID '_COMMON_FILTER_QC_ONLY.mat']), ...
  'FilterQC','ROIQC','TrialInventoryCompact', ...
  'cfgPilot','-v7.3');
end
writetable(FilterQC,fullfile(runDir,[runID '_FILTER_QC.csv']));
writetable(ROIQC,fullfile(runDir,[runID '_ROI_FILTER_QC.csv']));
writetable(Pilot.TrialInventory,fullfile(runDir,[runID '_TRIAL_INVENTORY.csv']));
fprintf('\nSource reconstruction complete: %s\n',runDir); disp(FilterQC); disp(ROIQC);
end

function D=crophalfopen(D,target)
for t=1:numel(D.trial)
 use=D.time{t}>=target(1) & D.time{t}<target(2);
 D.trial{t}=D.trial{t}(:,use); D.time{t}=D.time{t}(use);
 if isfield(D,'sampleinfo')&&size(D.sampleinfo,1)>=t
  first=find(use,1,'first'); last=find(use,1,'last'); old=D.sampleinfo(t,:);
  D.sampleinfo(t,:)=[old(1)+first-1 old(1)+last-1];
 end
end
end

function D=makeemgenvelope(D)
fs=double(D.fsample);
[b1,a1]=butter(4,[20 200]/(fs/2),'bandpass');
[b2,a2]=butter(4,10/(fs/2),'low');
for t=1:numel(D.trial)
 x=double(D.trial{t}); x=fillmissing(x,'linear',2,'EndValues','nearest');
 wide=filtfilt(b1,a1,x')'; rectified=abs(wide);
 D.trial{t}=filtfilt(b2,a2,rectified')';
end
end

function V=pickfield(S,names)
for k=1:numel(names), if isfield(S,names{k}), V=S.(names{k}); return; end, end
error('Expected variable missing: %s',strjoin(names,', '));
end

function F=getfilters(S)
if isfield(S,'avg')&&isfield(S.avg,'filter'), F=S.avg.filter;
elseif isfield(S,'filter'), F=S.filter;
else, error('LCMV output contains no filter field.'); end
end
