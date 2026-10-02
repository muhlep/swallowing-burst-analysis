function Results = burst_spatial_leakage_control(varargin)
% Leakage-aware beta ROI resolution analysis using the original forward
% models, common LCMV definition, and locked ROI-PCA weights.
%
% The input data and previous outputs are read-only. Subject-level results
% are saved after every participant, so a full-cohort run is resumable.

p = inputParser;
p.addParameter('sourceRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('originalCodeRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('fieldtripRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('outputRoot','',@(x)ischar(x)||isstring(x));
p.addParameter('resumeRunDir','',@(x)ischar(x)||isstring(x));
p.addParameter('subjectLimit',inf,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('lambda','5%',@(x)ischar(x)||isstring(x)||isnumeric(x));
p.addParameter('nPermutations',100000,@(x)isnumeric(x)&&isscalar(x)&&x>=9999);
p.addParameter('nBootstrap',20000,@(x)isnumeric(x)&&isscalar(x)&&x>=999);
p.addParameter('randomSeed',20260922,@(x)isnumeric(x)&&isscalar(x));
p.addParameter('runLabel','SPATIAL_LEAKAGE_CONTROL',@(x)ischar(x)||isstring(x));
p.addParameter('confirmRun',false,@(x)islogical(x)||ismember(x,[0 1]));
p.parse(varargin{:});
assert(p.Results.confirmRun,'Set confirmRun=true only after checking all paths.');

sourceRoot = char(string(p.Results.sourceRoot));
originalCodeRoot = char(string(p.Results.originalCodeRoot));
fieldtripRoot = char(string(p.Results.fieldtripRoot));
outputRoot = char(string(p.Results.outputRoot));
assert(isfolder(sourceRoot),'Source root missing: %s',sourceRoot);
assert(isfolder(originalCodeRoot),'Original code root missing: %s',originalCodeRoot);
assert(isfolder(fieldtripRoot),'FieldTrip root missing: %s',fieldtripRoot);
assert(isfolder(outputRoot),'Output root missing: %s',outputRoot);

addpath(originalCodeRoot);
addpath(fieldtripRoot);
ft_defaults;
assert(exist('pmd15_cfg','file')==2,'pmd15_cfg is unavailable.');

lockedHits = dir(fullfile(sourceRoot,'*_LOCKED_INPUTS.mat'));
assert(numel(lockedHits)==1,'Expected one top-level locked-input file.');
lockedFile = fullfile(lockedHits(1).folder,lockedHits(1).name);
L = load(lockedFile,'cfgFull','QCLocked','LockedSubjectQC', ...
    'EMGChannelMapping','TimeWindowDefinition');
subjectsAll = sort(string(L.LockedSubjectQC.Subject( ...
    ~L.LockedSubjectQC.ExcludeSubject)));
assert(numel(subjectsAll)==73,'Expected 73 healthy participants.');
nRequested = min(73,double(p.Results.subjectLimit));
subjects = subjectsAll(1:nRequested);

resumeRunDir = char(string(p.Results.resumeRunDir));
if isempty(resumeRunDir)
    stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
    runID = sprintf('BURST_%s_%s', ...
        char(string(p.Results.runLabel)),stamp);
    runDir = fullfile(outputRoot,runID);
    assert(~isfolder(runDir),'Refusing overwrite: %s',runDir);
    [ok,msg] = mkdir(runDir);
    assert(ok,'Cannot create run directory: %s',msg);
    cfg = struct('runID',runID,'runDir',runDir,'sourceRoot',sourceRoot, ...
        'originalCodeRoot',originalCodeRoot,'fieldtripRoot',fieldtripRoot, ...
        'subjects',subjects,'subjectLimit',nRequested,'band',[13 30], ...
        'bandLabel','beta','lambda',p.Results.lambda, ...
        'nPermutations',p.Results.nPermutations, ...
        'nBootstrap',p.Results.nBootstrap, ...
        'randomSeed',p.Results.randomSeed, ...
        'resolutionDefinition', ...
        'ROI virtual-sensor resolution matrix F*A using locked ROI-PC1 weights');
    Status = table(subjects,repmat("PENDING",nRequested,1), ...
        repmat("",nRequested,1),nan(nRequested,1), ...
        'VariableNames',{'Subject','State','Message','ElapsedSeconds'});
    save(fullfile(runDir,[runID '_STATE.mat']),'cfg','Status','-v7.3');
else
    runDir = resumeRunDir;
    assert(isfolder(runDir),'Resume directory missing: %s',runDir);
    stateHits = dir(fullfile(runDir,'*_STATE.mat'));
    assert(numel(stateHits)==1,'Expected one state file in resumeRunDir.');
    X = load(fullfile(stateHits(1).folder,stateHits(1).name),'cfg','Status');
    cfg = X.cfg;
    Status = X.Status;
    runID = cfg.runID;
    subjects = string(cfg.subjects);
    nRequested = numel(subjects);
    assert(strcmp(cfg.sourceRoot,sourceRoot),'Resume sourceRoot differs.');
end
stateFile = fullfile(runDir,[runID '_STATE.mat']);

betaFiles = dir(fullfile(sourceRoot,'**','*_ROI_PCA_PILOT.mat'));
allPaths = strings(numel(betaFiles),1);
for k = 1:numel(betaFiles)
    allPaths(k) = string(fullfile(betaFiles(k).folder,betaFiles(k).name));
end
betaFiles = betaFiles(contains(upper(allPaths),'_BETA_'));
assert(numel(betaFiles)==73,'Expected 73 beta ROI-PCA files.');

timerAll = tic;
for s = 1:nRequested
    sid = subjects(s);
    subjectFile = fullfile(runDir,sprintf('%s_SUBJECT_%s_RESOLUTION.mat', ...
        runID,char(sid)));
    if isfile(subjectFile)
        Status.State(s) = "COMPLETE";
        Status.Message(s) = "RESUMED_FROM_SUBJECT_FILE";
        fprintf('[%d/%d] %s already complete; skipped.\n',s,nRequested,sid);
        continue
    end
    fprintf('\n[%d/%d] BETA RESOLUTION: %s\n',s,nRequested,sid);
    subjectTimer = tic;
    try
        roiHit = betaFiles(contains( ...
            string({betaFiles.folder}) + filesep + string({betaFiles.name}), ...
            sid,'IgnoreCase',true));
        assert(numel(roiHit)==1,'Beta ROI file is not unique for %s.',sid);
        roiPath = fullfile(roiHit(1).folder,roiHit(1).name);
        R = load(roiPath,'ROIPilot','cfgROI');
        assert(isfield(R,'ROIPilot'),'ROIPilot missing for %s.',sid);
        assert(string(R.ROIPilot.bandLabel)=="beta", ...
            'Loaded ROI output is not beta for %s.',sid);

        [SubjectResolution,SubjectQC] = computeSubjectResolution( ...
            sid,L.QCLocked,L.LockedSubjectQC,L.EMGChannelMapping, ...
            L.TimeWindowDefinition,R.ROIPilot,p.Results.lambda);
        save(subjectFile,'SubjectResolution','SubjectQC','roiPath','-v7.3');
        writetable(SubjectResolution,fullfile(runDir,sprintf( ...
            '%s_SUBJECT_%s_EDGE_RESOLUTION.csv',runID,char(sid))));
        Status.State(s) = "COMPLETE";
        Status.Message(s) = "OK";
        Status.ElapsedSeconds(s) = toc(subjectTimer);
        save(stateFile,'cfg','Status','-v7.3');
        writetable(Status,fullfile(runDir,[runID '_STATUS.csv']));
        fprintf('%s complete | %.1f min | total %s\n',sid, ...
            Status.ElapsedSeconds(s)/60, ...
            char(duration(0,0,toc(timerAll),'Format','hh:mm:ss')));
    catch ME
        Status.State(s) = "FAILED";
        Status.Message(s) = string(ME.message);
        Status.ElapsedSeconds(s) = toc(subjectTimer);
        save(stateFile,'cfg','Status','ME','-v7.3');
        writetable(Status,fullfile(runDir,[runID '_STATUS.csv']));
        rethrow(ME)
    end
end

assert(all(Status.State=="COMPLETE"),'Not all requested participants completed.');
[SubjectEdgeResolution,SubjectQC] = collectSubjectFiles( ...
    runDir,runID,subjects);
GroupEdgeResolution = summarizeEdges(SubjectEdgeResolution, ...
    double(p.Results.nBootstrap),double(p.Results.randomSeed));
SetComparison = compareFinalEdgeSet(SubjectEdgeResolution, ...
    double(p.Results.nPermutations),double(p.Results.nBootstrap), ...
    double(p.Results.randomSeed));

Results = struct('InputQC',table(nRequested,nRequested==73, ...
    string(version),string(ft_version),string(lockedFile), ...
    'VariableNames',{'NSubjects','FullCohort','MATLABVersion', ...
    'FieldTripVersion','LockedInputFile'}), ...
    'SubjectQC',SubjectQC, ...
    'SubjectEdgeResolution',SubjectEdgeResolution, ...
    'GroupEdgeResolution',GroupEdgeResolution, ...
    'SetComparison',SetComparison, ...
    'cfg',cfg);

save(fullfile(runDir,[runID '_FINAL_RESULTS.mat']),'Results','-v7.3');
writetable(Results.InputQC,fullfile(runDir,[runID '_INPUT_QC.csv']));
writetable(SubjectQC,fullfile(runDir,[runID '_SUBJECT_QC.csv']));
writetable(SubjectEdgeResolution, ...
    fullfile(runDir,[runID '_SUBJECT_EDGE_RESOLUTION.csv']));
writetable(GroupEdgeResolution, ...
    fullfile(runDir,[runID '_GROUP_EDGE_RESOLUTION.csv']));
writetable(SetComparison, ...
    fullfile(runDir,[runID '_FINAL_EDGE_SET_COMPARISON.csv']));

fprintf('\nSpatial leakage input QC\n');
disp(Results.InputQC)
fprintf('\nSUBJECT QC\n');
disp(SubjectQC)
fprintf('\nLOCKED SIX BETA EDGES\n');
disp(GroupEdgeResolution(GroupEdgeResolution.LockedFinalBetaEdge,:))
fprintf('\nFINAL-EDGE SET COMPARISON\n');
disp(SetComparison)
fprintf('\nRESULT DIRECTORY\n%s\n',runDir);
end

function [E,Q] = computeSubjectResolution(sid,QCLocked,LockedSubjectQC, ...
    EMGChannelMapping,TimeWindowDefinition,ROIPilot,lambda)

conditions = ["EMG2";"EMG3"];
band = [13 30];
sr = string(LockedSubjectQC.Subject)==sid;
assert(sum(sr)==1 && ~LockedSubjectQC.ExcludeSubject(sr), ...
    'Participant is not uniquely included.');
mr = string(EMGChannelMapping.Subject)==sid;
assert(sum(mr)==1,'EMG mapping is not unique.');
dataset = char(EMGChannelMapping.Dataset(mr));
assert(isfolder(dataset),'Dataset missing: %s',dataset);

c = pmd15_cfg('subj',char(sid));
H = load(c.headmodel_file);
G = load(c.warped_grid_file);
R = load(c.roi_file,'ROIs');
vol = pickfield(H,{'vol','headmodel'});
grid = pickfield(G,{'warped_sourcemodel','sourcemodel'});
assert(isfield(R,'ROIs') && numel(R.ROIs)==21,'Expected 21 ROIs.');
vol = ft_convert_units(vol,'mm');
grid = ft_convert_units(grid,'mm');

TLavg = cell(2,1);
dataBand = cell(2,1);
nTrials = zeros(1,2);
for ci = 1:2
    cond = conditions(ci);
    wr = string(TimeWindowDefinition.Condition)==cond;
    assert(sum(wr)==1,'Window definition missing for %s.',cond);
    target = [TimeWindowDefinition.WindowStartSec(wr) ...
        TimeWindowDefinition.WindowEndSec(wr)];
    loadWindow = target + [-1 1];

    cfg0 = [];
    cfg0.dataset = dataset;
    cfg0.continuous = 'yes';
    cfg0.trialfun = 'ft_trialfun_general';
    cfg0.headerformat = 'ctf_ds';
    cfg0.dataformat = 'ctf_ds';
    cfg0.eventformat = 'ctf_ds';
    cfg0.trialdef.eventtype = char(cond);
    cfg0.trialdef.prestim = abs(loadWindow(1));
    cfg0.trialdef.poststim = loadWindow(2);
    cfg0 = ft_definetrial(cfg0);
    cfg0.channel = {'MEG'};
    cfg0.demean = 'no';
    D = ft_preprocessing(cfg0);
    nOriginal = numel(D.trial);
    D.trialinfo = (1:nOriginal)';

    qr = string(QCLocked.Subject)==sid & ...
        string(QCLocked.Condition)==cond;
    Qt = sortrows(QCLocked(qr,:),'OriginalTrial');
    assert(height(Qt)==nOriginal,'QC/raw trial mismatch for %s.',cond);
    assert(isequal(double(Qt.OriginalTrial),(1:nOriginal)'), ...
        'Original trial keys are not consecutive.');
    keep = ~Qt.ExcludeTrialGeneralArtifact;
    cfgs = [];
    cfgs.trials = find(keep);
    D = ft_selectdata(cfgs,D);

    cfgf = [];
    cfgf.channel = 'MEG';
    cfgf.bpfilter = 'yes';
    cfgf.bpfreq = band;
    cfgf.bpfilttype = 'but';
    cfgf.bpfiltdir = 'twopass';
    cfgf.demean = 'no';
    DB = ft_preprocessing(cfgf,D);
    DB = cropHalfOpen(DB,target);
    assert(all(cellfun(@(x)size(x,2)==600,DB.trial)), ...
        'Beta target trials must contain 600 samples.');
    cfgt = [];
    cfgt.covariance = 'yes';
    cfgt.covariancewindow = 'all';
    cfgt.keeptrials = 'no';
    TLavg{ci} = ft_timelockanalysis(cfgt,DB);
    dataBand{ci} = DB;
    nTrials(ci) = numel(DB.trial);
end

assert(isequal(TLavg{1}.label,TLavg{2}.label), ...
    'MEG channel labels differ across conditions.');
commonTL = TLavg{1};
commonTL.cov = (double(TLavg{1}.cov)+double(TLavg{2}.cov))/2;
commonTL.avg = (double(TLavg{1}.avg)+double(TLavg{2}.avg))/2;
if isfield(TLavg{1},'var') && isfield(TLavg{2},'var')
    commonTL.var = (double(TLavg{1}.var)+double(TLavg{2}.var))/2;
end

cfgsm = [];
cfgsm.sourcemodel = grid;
cfgsm.headmodel = vol;
cfgsm.grad = dataBand{1}.grad;
cfgsm.channel = cellstr(commonTL.label);
sourcemodel = ft_prepare_sourcemodel(cfgsm);
cfglf = [];
cfglf.headmodel = vol;
cfglf.sourcemodel = sourcemodel;
cfglf.grad = dataBand{1}.grad;
cfglf.channel = cellstr(commonTL.label);
cfglf.reducerank = 2;
cfglf.normalize = 'yes';
leadfield = ft_prepare_leadfield(cfglf,commonTL);

cfgsa = [];
cfgsa.method = 'lcmv';
cfgsa.sourcemodel = leadfield;
cfgsa.headmodel = vol;
cfgsa.lcmv.keepfilter = 'yes';
cfgsa.lcmv.fixedori = 'yes';
cfgsa.lcmv.projectmom = 'yes';
cfgsa.lcmv.lambda = lambda;
sourceCommon = ft_sourceanalysis(cfgsa,commonTL);
filters = getFilters(sourceCommon);

roiLabels = string(ROIPilot.roiLabels(:));
assert(numel(roiLabels)==21,'Locked ROI-PCA output must contain 21 ROIs.');
assert(numel(ROIPilot.roiWeights)==21 && ...
    numel(ROIPilot.roiIndices)==21,'Locked ROI weights/indices missing.');
[F,A] = makeROIOperators(filters,leadfield,ROIPilot);
resolution = F*A;
diagonalGain = diag(resolution);
assert(all(isfinite(diagonalGain)) && all(abs(diagonalGain)>1e-12), ...
    'ROI resolution matrix has invalid diagonal gain.');

sourceNormalized = abs(resolution ./ diagonalGain');
targetNormalized = abs(resolution ./ diagonalGain);
reciprocal = sqrt(sourceNormalized .* sourceNormalized');
forwardSimilarity = abs(cosineColumns(A));
filterSimilarity = abs(cosineRows(F));

[ia,ib] = find(triu(true(21),1));
nEdge = numel(ia);
rows = cell(nEdge,12);
for e = 1:nEdge
    a = ia(e); b = ib(e);
    rows(e,:) = {sid,roiLabels(a),roiLabels(b),a,b, ...
        sourceNormalized(a,b),sourceNormalized(b,a), ...
        max(sourceNormalized(a,b),sourceNormalized(b,a)), ...
        reciprocal(a,b),forwardSimilarity(a,b),filterSimilarity(a,b), ...
        isLockedFinalEdge(roiLabels(a),roiLabels(b))};
end
E = cell2table(rows,'VariableNames',{'Subject','ROI_A','ROI_B', ...
    'ROI_A_Index','ROI_B_Index','Leakage_B_to_A_SourceNormalized', ...
    'Leakage_A_to_B_SourceNormalized','MaximumDirectionalLeakage', ...
    'ReciprocalGeometricLeakage','ForwardTopographySimilarity', ...
    'InverseFilterSimilarity','LockedFinalBetaEdge'});

off = ~eye(21);
Q = table(sid,nTrials(1),nTrials(2),rcond(double(commonTL.cov)), ...
    min(abs(diagonalGain)),max(abs(diagonalGain)), ...
    median(reciprocal(off),'omitnan'),max(reciprocal(off)), ...
    median(forwardSimilarity(off),'omitnan'),max(forwardSimilarity(off)), ...
    'VariableNames',{'Subject','NTrialsEMG2','NTrialsEMG3', ...
    'CovarianceRcond','MinimumAbsoluteDiagonalGain', ...
    'MaximumAbsoluteDiagonalGain','MedianOffDiagonalReciprocalLeakage', ...
    'MaximumOffDiagonalReciprocalLeakage', ...
    'MedianOffDiagonalForwardSimilarity','MaximumOffDiagonalForwardSimilarity'});
end

function [F,A] = makeROIOperators(filters,leadfield,ROIPilot)
nROI = 21;
nChannel = numel(filters{find(~cellfun(@isempty,filters),1)});
F = zeros(nROI,nChannel);
A = zeros(nChannel,nROI);
for r = 1:nROI
    idx = double(ROIPilot.roiIndices{r}(:));
    weights = double(ROIPilot.roiWeights{r}(:));
    assert(numel(idx)==numel(weights),'ROI index/weight mismatch.');
    for v = 1:numel(idx)
        j = idx(v);
        assert(j>=1 && j<=numel(filters) && ...
            ~isempty(filters{j}),'Missing spatial filter in ROI.');
        assert(j<=numel(leadfield.leadfield) && ...
            ~isempty(leadfield.leadfield{j}),'Missing leadfield in ROI.');
        wj = double(filters{j}(:)');
        Lj = double(leadfield.leadfield{j});
        assert(size(Lj,1)==nChannel && size(wj,2)==nChannel, ...
            'Leadfield/filter channel mismatch.');
        pass = wj*Lj;
        assert(norm(pass)>0,'Cannot determine fixed orientation.');
        q = pass'/norm(pass);
        topography = Lj*q;
        F(r,:) = F(r,:) + weights(v)*wj;
        A(:,r) = A(:,r) + weights(v)*topography;
    end
end
assert(all(isfinite(F),'all') && all(isfinite(A),'all'), ...
    'Non-finite ROI operators.');
end

function [E,Q] = collectSubjectFiles(runDir,runID,subjects)
edgeParts = cell(numel(subjects),1);
qcParts = cell(numel(subjects),1);
for s = 1:numel(subjects)
    f = fullfile(runDir,sprintf('%s_SUBJECT_%s_RESOLUTION.mat', ...
        runID,char(subjects(s))));
    assert(isfile(f),'Subject result missing: %s',f);
    X = load(f,'SubjectResolution','SubjectQC');
    edgeParts{s} = X.SubjectResolution;
    qcParts{s} = X.SubjectQC;
end
E = vertcat(edgeParts{:});
Q = vertcat(qcParts{:});
assert(height(E)==210*numel(subjects),'Unexpected subject-edge row count.');
end

function G = summarizeEdges(E,nBootstrap,seed)
edgeKeys = unique(E(:,{'ROI_A','ROI_B'}),'rows','stable');
nEdge = height(edgeKeys);
rows = cell(nEdge,16);
rng(seed,'twister');
for e = 1:nEdge
    q = string(E.ROI_A)==string(edgeKeys.ROI_A(e)) & ...
        string(E.ROI_B)==string(edgeKeys.ROI_B(e));
    X = E(q,:);
    x = double(X.ReciprocalGeometricLeakage);
    y = double(X.MaximumDirectionalLeakage);
    z = double(X.ForwardTopographySimilarity);
    ciX = bootstrapMeanCI(x,nBootstrap);
    ciY = bootstrapMeanCI(y,nBootstrap);
    ciZ = bootstrapMeanCI(z,nBootstrap);
    rows(e,:) = {string(edgeKeys.ROI_A(e)),string(edgeKeys.ROI_B(e)), ...
        height(X),mean(x),median(x),ciX(1),ciX(2), ...
        mean(y),median(y),ciY(1),ciY(2), ...
        mean(z),median(z),ciZ(1),ciZ(2),logical(X.LockedFinalBetaEdge(1))};
end
G = cell2table(rows,'VariableNames',{'ROI_A','ROI_B','NSubjects', ...
    'MeanReciprocalGeometricLeakage','MedianReciprocalGeometricLeakage', ...
    'ReciprocalMeanCI95Lower','ReciprocalMeanCI95Upper', ...
    'MeanMaximumDirectionalLeakage','MedianMaximumDirectionalLeakage', ...
    'MaximumDirectionalMeanCI95Lower','MaximumDirectionalMeanCI95Upper', ...
    'MeanForwardTopographySimilarity','MedianForwardTopographySimilarity', ...
    'ForwardSimilarityMeanCI95Lower','ForwardSimilarityMeanCI95Upper', ...
    'LockedFinalBetaEdge'});
G.ReciprocalLeakagePercentileAmong210 = tiedPercentile( ...
    double(G.MeanReciprocalGeometricLeakage));
G.MaximumDirectionalPercentileAmong210 = tiedPercentile( ...
    double(G.MeanMaximumDirectionalLeakage));
G.ForwardSimilarityPercentileAmong210 = tiedPercentile( ...
    double(G.MeanForwardTopographySimilarity));
end

function T = compareFinalEdgeSet(E,nPerm,nBoot,seed)
subjects = unique(string(E.Subject),'stable');
n = numel(subjects);
d = nan(n,3);
metrics = {'ReciprocalGeometricLeakage', ...
    'MaximumDirectionalLeakage','ForwardTopographySimilarity'};
for s = 1:n
    X = E(string(E.Subject)==subjects(s),:);
    final = logical(X.LockedFinalBetaEdge);
    assert(sum(final)==6,'Expected six locked beta edges.');
    for m = 1:3
        values = double(X.(metrics{m}));
        d(s,m) = mean(values(final))-mean(values(~final));
    end
end
rng(seed,'twister');
rows = cell(3,10);
for m = 1:3
    x = d(:,m);
    mu = mean(x);
    ci = bootstrapMeanCI(x,nBoot);
    dz = mu/std(x,0);
    if n==73
        pValue = signFlipP(x,nPerm);
    else
        pValue = NaN;
    end
    rows(m,:) = {string(metrics{m}),n,meanSet(E,metrics{m},true), ...
        meanSet(E,metrics{m},false),mu,ci(1),ci(2),dz,pValue,n==73};
end
T = cell2table(rows,'VariableNames',{'Metric','NSubjects', ...
    'MeanLockedSix','MeanOther204','LockedMinusOther', ...
    'BootstrapCI95Lower','BootstrapCI95Upper','CohenDz', ...
    'TwoSidedSignFlipP','ConfirmatoryFullCohort'});
end

function value = meanSet(E,metric,finalFlag)
value = mean(double(E.(metric)(logical(E.LockedFinalBetaEdge)==finalFlag)));
end

function ci = bootstrapMeanCI(x,nBoot)
n = numel(x);
boot = zeros(nBoot,1);
done = 0;
while done<nBoot
    m = min(2000,nBoot-done);
    idx = randi(n,n,m);
    sample = reshape(x(idx),n,m);
    boot(done+(1:m)) = mean(sample,1)';
    done = done+m;
end
ci = prctile(boot,[2.5 97.5]);
end

function pValue = signFlipP(x,nPerm)
n = numel(x);
observed = abs(mean(x));
count = 0;
done = 0;
while done<nPerm
    m = min(5000,nPerm-done);
    signs = 2*(rand(m,n)>0.5)-1;
    nullMeans = abs((signs*x)/n);
    count = count + sum(nullMeans>=observed);
    done = done+m;
end
pValue = (count+1)/(nPerm+1);
end

function p = tiedPercentile(x)
n = numel(x);
p = nan(n,1);
for k = 1:n
    p(k) = 100*(sum(x<x(k)) + 0.5*sum(x==x(k)))/n;
end
end

function tf = isLockedFinalEdge(a,b)
pairs = [ ...
    "BA44_R","Ins_Ant_R"; ...
    "Ins_Ant_L","Ins_Post_L"; ...
    "Ins_Ant_R","Ins_Post_R"; ...
    "Ins_Post_R","PMd_R"; ...
    "Ins_Post_R","S2_R"; ...
    "Ins_Post_R","SMG_R"];
tf = false;
for k = 1:size(pairs,1)
    tf = tf || ((a==pairs(k,1) && b==pairs(k,2)) || ...
        (a==pairs(k,2) && b==pairs(k,1)));
end
end

function C = cosineColumns(X)
n = size(X,2);
C = eye(n);
for a = 1:n
    for b = a+1:n
        value = dot(X(:,a),X(:,b))/(norm(X(:,a))*norm(X(:,b)));
        C(a,b) = value;
        C(b,a) = value;
    end
end
end

function C = cosineRows(X)
C = cosineColumns(X');
end

function D = cropHalfOpen(D,target)
for t = 1:numel(D.trial)
    use = D.time{t}>=target(1) & D.time{t}<target(2);
    D.trial{t} = D.trial{t}(:,use);
    D.time{t} = D.time{t}(use);
    if isfield(D,'sampleinfo') && size(D.sampleinfo,1)>=t
        first = find(use,1,'first');
        last = find(use,1,'last');
        old = D.sampleinfo(t,:);
        D.sampleinfo(t,:) = [old(1)+first-1 old(1)+last-1];
    end
end
end

function V = pickfield(S,names)
for k = 1:numel(names)
    if isfield(S,names{k})
        V = S.(names{k});
        return
    end
end
error('Expected variable missing: %s',strjoin(names,', '));
end

function F = getFilters(S)
if isfield(S,'avg') && isfield(S.avg,'filter')
    F = S.avg.filter;
elseif isfield(S,'filter')
    F = S.filter;
else
    error('LCMV output contains no filter field.');
end
end
