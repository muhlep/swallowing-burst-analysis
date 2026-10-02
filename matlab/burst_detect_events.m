function [Events,Inventory,Normalization] = burst_detect_events( ...
    X,conditions,trialKeys,labels,bandLabel,fs,varargin)
% Detect nested threshold events in two-condition time-series data.
% X is a 2-element cell; each element is signal x trial x sample.
% Normalization is learned jointly from EMG2 and EMG3 for every signal.
% Higher-threshold events are children of z=1.5 parent excursions. This
% propagates left censoring consistently and guarantees nested prevalence.

p=inputParser;
p.addParameter('thresholds',[1.5 2.0 2.5],@(x)isnumeric(x)&&isvector(x));
p.addParameter('minimumDurationSamples',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('signalType','CORTICAL',@(x)ischar(x)||isstring(x));
p.parse(varargin{:});

thresholds=sort(double(p.Results.thresholds(:)'));
assert(isequal(thresholds,[1.5 2 2.5]),'Required thresholds are [1.5 2 2.5].');
minSamples=ceil(double(p.Results.minimumDurationSamples));
conditions=string(conditions(:)); labels=string(labels(:));
assert(numel(X)==2&&numel(conditions)==2&&numel(trialKeys)==2,'Two conditions required.');
assert(size(X{1},1)==numel(labels)&&size(X{2},1)==numel(labels),'Label count mismatch.');
assert(size(X{1},3)==600&&size(X{2},3)==600,'Expected 600 target samples.');

nSignal=numel(labels); nNorm=nSignal;
normRows=cell(nNorm,8); eventRows=cell(0,19); inventoryRows=cell(0,16);
eventID=0;

for r=1:nSignal
 pooled=[];
 for ciNorm=1:2
  for tNorm=1:size(X{ciNorm},2)
   xNorm=squeeze(double(X{ciNorm}(r,tNorm,:)));
   if string(p.Results.signalType)=="CORTICAL"
    xNorm=abs(hilbert(xNorm));
   end
   pooled=[pooled;xNorm(:)]; %#ok<AGROW>
  end
 end
 center=median(pooled,'omitnan');
 scale=1.4826*median(abs(pooled-center),'omitnan');
 fallback=false;
 if ~isfinite(scale)||scale<=eps(max(1,abs(center)))
  scale=std(pooled,0,'omitnan'); fallback=true;
 end
 assert(isfinite(center)&&isfinite(scale)&&scale>0,'Invalid robust normalization for %s.',labels(r));
 normRows(r,:)={labels(r),string(bandLabel),center,scale,numel(pooled), ...
  sum(isfinite(pooled)),fallback,string(p.Results.signalType)};

 for ci=1:2
  nTrial=size(X{ci},2); keys=double(trialKeys{ci}(:));
  assert(numel(keys)==nTrial,'Trial-key count mismatch for %s.',conditions(ci));
  for t=1:nTrial
   raw=squeeze(double(X{ci}(r,t,:)))';
   assert(numel(raw)==600&&all(isfinite(raw)),'Invalid trial signal.');
   if string(p.Results.signalType)=="CORTICAL"
    amplitude=abs(hilbert(raw));
   else
    amplitude=raw;
   end
   z=(amplitude-center)./scale;

   % Parent excursions use the lowest threshold before duration filtering.
   parents=findSegments(z>=thresholds(1));
   parentId=zeros(1,numel(z)); parentLeft=false(size(parents,1),1);
   for q=1:size(parents,1)
    parentId(parents(q,1):parents(q,2))=q;
    parentLeft(q)=parents(q,1)==1;
   end

   for h=1:numel(thresholds)
    thr=thresholds(h); seg=findSegments(z>=thr);
    keep=(seg(:,2)-seg(:,1)+1)>=minSamples; seg=seg(keep,:);
    nAll=size(seg,1); nPrimary=0; anyLeft=false; anyRight=false;
    for e=1:nAll
     a=seg(e,1); b=seg(e,2); pid=parentId(a);
     assert(pid>0,'Threshold event lacks z=1.5 parent.');
     left=parentLeft(pid); right=(b==numel(z)); primary=~left;
     nPrimary=nPrimary+double(primary); anyLeft=anyLeft||left; anyRight=anyRight||right;
     zz=z(a:b); [peakZ,peakRelative]=max(zz); peakSample=a+peakRelative-1;
     eventID=eventID+1;
     eventRows(end+1,:)={eventID,string(p.Results.signalType),labels(r),string(bandLabel), ...
      conditions(ci),keys(t),t,thr,pid,a,b,peakSample,(b-a+1), ...
      (b-a+1)/fs,peakZ,sum(max(zz-thr,0))/fs,left,right,primary}; %#ok<AGROW>
    end
    inventoryRows(end+1,:)={string(p.Results.signalType),labels(r),string(bandLabel), ...
     conditions(ci),keys(t),t,thr,true,nAll>0,nPrimary>0,nAll,nPrimary, ...
     anyLeft,anyRight,minSamples,minSamples/fs}; %#ok<AGROW>
   end
  end
 end
end

Events=cell2table(eventRows,'VariableNames',{'EventID','SignalType','Signal','Band', ...
 'Condition','OriginalTrial','TrialIndex','ZThreshold','ParentID','OnsetSample', ...
 'OffsetSample','PeakSample','DurationSamples','DurationSec','PeakZ','AUCAboveThreshold', ...
 'LeftCensored','RightCensored','PrimaryEligible'});
Inventory=cell2table(inventoryRows,'VariableNames',{'SignalType','Signal','Band', ...
 'Condition','OriginalTrial','TrialIndex','ZThreshold','TrialValid','AnyBurstInclusive', ...
 'AnyBurstPrimary','NBurstsInclusive','NBurstsPrimary','AnyLeftCensored', ...
 'AnyRightCensored','MinimumDurationSamples','MinimumDurationSec'});
Normalization=cell2table(normRows,'VariableNames',{'Signal','Band','NormCenter', ...
 'NormScale','NValues','NFinite','ScaleFallback','SignalType'});
end

function seg=findSegments(mask)
mask=logical(mask(:)'); d=diff([false mask false]);
seg=[find(d==1)' find(d==-1)'-1];
end
