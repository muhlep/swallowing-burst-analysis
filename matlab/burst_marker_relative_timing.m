function Results = burst_marker_relative_timing(burstMatFile,varargin)
% One participant-synchronous max-|T| family across four window widths,
% three standard bands and three direct band contrasts at z=2.

p=inputParser;
addParameter(p,'variableName','BurstFullR02',@(x)ischar(x)||isstring(x));
addParameter(p,'outputRoot','',@(x)ischar(x)||isstring(x));
addParameter(p,'nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
addParameter(p,'nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
addParameter(p,'seed',20260923,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'confirmRun',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
assert(p.Results.confirmRun,'Set confirmRun=true after checking the input file.');

X=load(char(string(burstMatFile)),char(string(p.Results.variableName)));
B=X.(char(string(p.Results.variableName))); E=B.CorticalEvents; I=B.CorticalInventory;
subjects=sort(unique(string(I.Subject))); bands=["theta";"alpha";"beta"];
widths=[0.1 0.2 0.3 0.4]; fs=600;
assert(numel(subjects)==73,'Expected 73 participants.');
D=nan(73,12); Pre=nan(73,12); Post=nan(73,12); endpoint=0;
bandLabel=strings(12,1); widthLabel=nan(12,1);
for w=1:4
    for b=1:3
        endpoint=endpoint+1; band=bands(b); width=widths(w);
        bandLabel(endpoint)=band; widthLabel(endpoint)=width;
        iq=string(I.Condition)=="EMG2" & string(I.Band)==band & ...
            double(I.ZThreshold)==2 & logical(I.TrialValid);
        eq=string(E.Condition)=="EMG2" & string(E.Band)==band & ...
            double(E.ZThreshold)==2 & logical(E.PrimaryEligible);
        IT=I(iq,:); ET=E(eq,:);
        keyI=string(IT.Subject)+"|"+string(IT.Signal)+"|"+string(IT.OriginalTrial);
        keyE=string(ET.Subject)+"|"+string(ET.Signal)+"|"+string(ET.OriginalTrial);
        [ok,loc]=ismember(keyE,keyI); assert(all(ok),'Event/inventory mismatch.');
        onset=-0.4+(double(ET.OnsetSample)-1)/fs;
        pre=false(height(IT),1); post=false(height(IT),1);
        for k=1:height(ET)
            pre(loc(k))=pre(loc(k)) || (onset(k)>=-width && onset(k)<0);
            post(loc(k))=post(loc(k)) || (onset(k)>=0 && onset(k)<width);
        end
        for s=1:73
            q=string(IT.Subject)==subjects(s);
            Pre(s,endpoint)=100*mean(pre(q)); Post(s,endpoint)=100*mean(post(q));
            D(s,endpoint)=Post(s,endpoint)-Pre(s,endpoint);
        end
    end
end

% Add alpha-theta, beta-theta and beta-alpha for every width.
C=nan(73,12); contrastLabel=strings(12,1); contrastWidth=nan(12,1); e=0;
for w=1:4
    idx=(w-1)*3+(1:3); Xw=D(:,idx);
    names=["alpha - theta";"beta - theta";"beta - alpha"];
    weights=[-1 1 0;-1 0 1;0 -1 1];
    for c=1:3
        e=e+1; C(:,e)=Xw*weights(c,:)'; contrastLabel(e)=names(c); contrastWidth(e)=widths(w);
    end
end
All=[D C];
rng(double(p.Results.seed),'twister');
[rawP,maxTP,tObs]=jointSignFlip(All,double(p.Results.nPermutations));

rows=cell(24,14);
for j=1:24
    x=All(:,j); n=numel(x); boot=zeros(double(p.Results.nBootstrap),1); done=0;
    while done<numel(boot)
        m=min(2000,numel(boot)-done); idx=randi(n,n,m);
        boot(done+(1:m))=mean(reshape(x(idx),n,m),1)'; done=done+m;
    end
    ci=prctile(boot,[2.5 97.5]);
    if j<=12
        type="within-band post minus pre"; label=bandLabel(j); width=widthLabel(j);
        preMean=mean(Pre(:,j)); postMean=mean(Post(:,j));
    else
        type="difference of band-specific post-minus-pre effects";
        label=contrastLabel(j-12); width=contrastWidth(j-12); preMean=NaN; postMean=NaN;
    end
    rows(j,:)={type,label,width,73,preMean,postMean,mean(x),ci(1),ci(2), ...
        mean(x)/std(x,0),tObs(j),rawP(j),maxTP(j),mean(x>0)};
end
Stats=cell2table(rows,'VariableNames',{'EndpointType','BandOrContrast','HalfWidthSec', ...
    'NSubjects','MeanPrePercent','MeanPostPercent','EstimatePP', ...
    'BootstrapCI95Lower','BootstrapCI95Upper','CohenDz','ObservedT','RawP', ...
    'Unified24EndpointMaxTFWERP','ProportionPositive'});

subjectRows=cell(73*24,6); row=0;
for j=1:24
    for s=1:73
        row=row+1; subjectRows(row,:)={subjects(s),rows{j,1},rows{j,2},rows{j,3},All(s,j),2};
    end
end
SubjectData=cell2table(subjectRows,'VariableNames',{'Subject','EndpointType', ...
    'BandOrContrast','HalfWidthSec','EstimatePP','ZThreshold'});
runID="BURST_MARKER_RELATIVE_TIMING_"+ ...
    string(datetime('now','Format','yyyyMMdd_HHmmss'));
runDir=fullfile(char(string(p.Results.outputRoot)),char(runID));
assert(~isfolder(runDir),'Refusing overwrite.'); mkdir(runDir);
writetable(Stats,fullfile(runDir,'UNIFIED_TIMING_24_ENDPOINT_STATS.csv'));
writetable(SubjectData,fullfile(runDir,'UNIFIED_TIMING_PARTICIPANT_DATA.csv'));
cfg=struct('runID',runID,'runDir',string(runDir),'inputFile',string(burstMatFile), ...
    'family','4 widths x (3 bands + 3 direct band contrasts) at z=2', ...
    'nEndpoints',24,'nPermutations',p.Results.nPermutations, ...
    'nBootstrap',p.Results.nBootstrap,'seed',p.Results.seed,'status','EXPLORATORY');
Results=struct('Stats',Stats,'SubjectData',SubjectData,'cfg',cfg);
save(fullfile(runDir,char(runID+"_RESULTS.mat")),'Results','cfg','-v7.3');
fprintf('\nMarker-relative timing analysis complete\n%s\n',runDir); disp(Stats)
end

function [rawP,maxTP,tObs]=jointSignFlip(D,nPerm)
n=size(D,1); tObs=mean(D,1)./(std(D,0,1)/sqrt(n));
rawCount=zeros(1,size(D,2)); maxCount=zeros(1,size(D,2)); done=0;
while done<nPerm
    m=min(2000,nPerm-done); signs=2*(rand(n,m)>.5)-1;
    mu=(signs'*D)/n; sd=nan(m,size(D,2));
    for j=1:size(D,2), sd(:,j)=std(signs.*D(:,j),0,1)'; end
    nullT=mu./(sd/sqrt(n)); mx=max(abs(nullT),[],2);
    rawCount=rawCount+sum(abs(nullT)>=abs(tObs),1);
    for j=1:size(D,2), maxCount(j)=maxCount(j)+sum(mx>=abs(tObs(j))); end
    done=done+m;
end
rawP=(rawCount+1)/(nPerm+1); maxTP=(maxCount+1)/(nPerm+1);
end
