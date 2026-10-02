function Results = burst_shared_roi_morphology(burstMatFile,varargin)
% Shared-ROI sensitivity for duration, peak and suprathreshold area.
% Within each participant and band, only ROIs with complete-event estimates
% in both conditions contribute to either condition summary.

p=inputParser;
addParameter(p,'variableName','BurstFullR02',@(x)ischar(x)||isstring(x));
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
addParameter(p,'seed',20260923,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun,'Set confirmRun=true after checking the input file.');

burstMatFile=char(string(burstMatFile));
assert(isfile(burstMatFile),'Input MAT file not found.');
X=load(burstMatFile,char(string(p.Results.variableName)));
B=X.(char(string(p.Results.variableName)));
E=B.CorticalEvents;
required={'Subject','Condition','Band','Signal','ZThreshold','PrimaryEligible', ...
    'RightCensored','DurationSec','PeakZ','AUCAboveThreshold'};
assert(all(ismember(required,E.Properties.VariableNames)),'Event schema mismatch.');

subjects=sort(unique(string(E.Subject)));
bands=["theta";"alpha";"beta"];
conditions=["EMG2";"EMG3"];
metrics=["DurationSec";"PeakZ";"AUCAboveThreshold"];
assert(numel(subjects)==73,'Expected 73 participants.');

rows=cell(73*3*3,13); row=0;
D=nan(73,9); endpointBand=strings(9,1); endpointMetric=strings(9,1); endpoint=0;
for b=1:numel(bands)
    band=bands(b);
    bandMask=string(E.Band)==band;
    roiNames=sort(unique(string(E.Signal(bandMask))));
    assert(numel(roiNames)==21,'Expected 21 ROIs for %s.',band);
    V=nan(73,21,2,3);
    for s=1:73
        for r=1:21
            for c=1:2
                q=string(E.Subject)==subjects(s) & string(E.Band)==band & ...
                    string(E.Signal)==roiNames(r) & string(E.Condition)==conditions(c) & ...
                    double(E.ZThreshold)==2 & logical(E.PrimaryEligible) & ...
                    ~logical(E.RightCensored);
                for m=1:3
                    metricName=char(metrics(m));
                    metricValues=E.(metricName);
                    x=double(metricValues(q));
                    if ~isempty(x), V(s,r,c,m)=median(x,'omitnan'); end
                end
            end
        end
    end
    for m=1:3
        endpoint=endpoint+1;
        endpointBand(endpoint)=band; endpointMetric(endpoint)=metrics(m);
        for s=1:73
            finite1=isfinite(V(s,:,1,m)); finite2=isfinite(V(s,:,2,m));
            shared=finite1 & finite2;
            emg2=NaN; emg3=NaN; difference=NaN;
            if any(shared)
                emg2=median(reshape(V(s,shared,1,m),[],1),'omitnan');
                emg3=median(reshape(V(s,shared,2,m),[],1),'omitnan');
                difference=emg2-emg3;
            end
            D(s,endpoint)=difference;
            row=row+1;
            rows(row,:)={subjects(s),band,metrics(m),sum(finite1),sum(finite2), ...
                sum(shared),emg2,emg3,difference,2, ...
                "complete events only","median within ROI then across shared ROIs", ...
                "EMG2 minus EMG3"};
        end
    end
end
SubjectStats=cell2table(rows,'VariableNames',{'Subject','Band','Metric', ...
    'NContributingROIsEMG2','NContributingROIsEMG3','NSharedROIs','EMG2', ...
    'EMG3','Difference','ZThreshold','EventSupport','Aggregation','Contrast'});

rng(double(p.Results.seed),'twister');
[rawP,maxTP,tObs]=jointSignFlipFinite(D,double(p.Results.nPermutations));
summaryRows=cell(9,13);
for e=1:9
    d=D(:,e); keep=isfinite(d); x=d(keep); n=numel(x);
    boot=zeros(double(p.Results.nBootstrap),1); done=0;
    while done<numel(boot)
        q=min(2000,numel(boot)-done); idx=randi(n,n,q);
        boot(done+(1:q))=mean(reshape(x(idx),n,q),1)'; done=done+q;
    end
    ci=prctile(boot,[2.5 97.5]);
    sub=SubjectStats(string(SubjectStats.Band)==endpointBand(e) & ...
        string(SubjectStats.Metric)==endpointMetric(e),:);
    summaryRows(e,:)={endpointBand(e),endpointMetric(e),n, ...
        mean(sub.EMG2,'omitnan'),mean(sub.EMG3,'omitnan'),mean(x),ci(1),ci(2), ...
        mean(x)/std(x,0),tObs(e),rawP(e),maxTP(e),median(sub.NSharedROIs,'omitnan')};
end
Summary=cell2table(summaryRows,'VariableNames',{'Band','Metric','NPaired', ...
    'MeanEMG2','MeanEMG3','MeanDifference','BootstrapCI95Lower', ...
    'BootstrapCI95Upper','CohenDz','ObservedT','RawP','JointMaxTFWERP', ...
    'MedianSharedROIs'});

runID="BURST_SHARED_ROI_MORPHOLOGY_"+ ...
    string(datetime('now','Format','yyyyMMdd_HHmmss'));
runDir=fullfile(char(string(p.Results.outputRoot)),char(runID));
assert(~isfolder(runDir),'Refusing overwrite.'); mkdir(runDir);
writetable(SubjectStats,fullfile(runDir,'PARTICIPANT_SHARED_ROI_MORPHOLOGY.csv'));
writetable(Summary,fullfile(runDir,'SHARED_ROI_MORPHOLOGY_STATS.csv'));
cfg=struct('runID',runID,'runDir',string(runDir),'inputFile',string(burstMatFile), ...
    'zThreshold',2,'family','3 bands x 3 morphology endpoints', ...
    'nPermutations',p.Results.nPermutations,'nBootstrap',p.Results.nBootstrap, ...
    'seed',p.Results.seed,'status','SENSITIVITY');
Results=struct('SubjectStats',SubjectStats,'Summary',Summary,'cfg',cfg);
save(fullfile(runDir,char(runID+"_RESULTS.mat")),'Results','cfg','-v7.3');
fprintf('\nShared-ROI morphology analysis complete\n%s\n',runDir); disp(Summary)
end

function [rawP,maxTP,tObs]=jointSignFlipFinite(D,nPerm)
n=size(D,1); p=size(D,2); tObs=nan(1,p);
for j=1:p, x=D(:,j); x=x(isfinite(x)); tObs(j)=mean(x)/(std(x,0)/sqrt(numel(x))); end
rawCount=zeros(1,p); maxCount=zeros(1,p); done=0;
while done<nPerm
    m=min(2000,nPerm-done); signs=2*(rand(n,m)>.5)-1; nullT=nan(m,p);
    for j=1:p
        keep=isfinite(D(:,j)); x=D(keep,j); S=signs(keep,:);
        nullT(:,j)=(mean(S.*x,1)./(std(S.*x,0,1)/sqrt(numel(x))))';
        rawCount(j)=rawCount(j)+sum(abs(nullT(:,j))>=abs(tObs(j)));
    end
    mx=max(abs(nullT),[],2);
    for j=1:p, maxCount(j)=maxCount(j)+sum(mx>=abs(tObs(j))); end
    done=done+m;
end
rawP=(rawCount+1)/(nPerm+1); maxTP=(maxCount+1)/(nPerm+1);
end
