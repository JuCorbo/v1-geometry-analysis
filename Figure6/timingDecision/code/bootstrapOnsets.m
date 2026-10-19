%% Bootstrap CIs on the pooled naive-vs-trained divergence ONSET, per geometric measure
% Resample sessions (with replacement) within each group, rerun groupSeparationTest,
% record onset -> percentile CI. Justified pooling (day-to-day dynamics similar).
clear; clc;
codeDir = fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
anaDir = fileparts(codeDir);
outDir = codeDir;
figRoot = fileparts(fileparts(anaDir));
addpath(fullfile(figRoot,'_sharedWindowedData','code'));  % groupSeparationTest
rng(7);

Nboot  = 1000;
nPerms = 1000;
onsetOffsetSec = 1*0.0646;  % 1*dt: report at window centre (tvec_centre)

%% ---- assemble per-session pooled data for each LOCAL measure ----
% dimensionality (D_PR)
D = load(fullfile(anaDir,'cache','dimTiming.mat'));           % dimN,dimT {1x6}, tvec
tvec = D.tvec(:)';
meas(1).name='dimDPR';   meas(1).N=vertcat(D.dimN{:}); meas(1).T=vertcat(D.dimT{:}); meas(1).tail='left';  meas(1).minSust=1;
% radius & centroid distance
S = load(fullfile(figRoot,'Figure3','nativeSpaceGeometry','cache','originalSpaceBasicGeometrySession.mat'),'radius_bySess','centDist_bySess','tvec');  % 0.12 uniformized geometry
poolGrp = @(cb,row) cell2mat(cellfun(@(c) c(:)', [cb{row,:}]', 'uni',0));
meas(2).name='radiusMix'; meas(2).N=poolGrp(S.radius_bySess,1);   meas(2).T=poolGrp(S.radius_bySess,2);   meas(2).tail='left';  meas(2).minSust=1;
meas(3).name='centDist';  meas(3).N=poolGrp(S.centDist_bySess,1); meas(3).T=poolGrp(S.centDist_bySess,2); meas(3).tail='right'; meas(3).minSust=1;

assert(numel(S.tvec)==numel(tvec) && all(abs(S.tvec(:)'-tvec)<1e-6), 'tvec mismatch');

%% ---- point onset + bootstrap ----
results = struct();
for m = 1:numel(meas)
    N = meas(m).N; T = meas(m).T; tail = meas(m).tail; ms = meas(m).minSust;
    nN=size(N,1); nT=size(T,1);
    gargs = {'tail',tail,'correction','cluster','alpha',0.05,'minSustained',ms, ...
             'minHedgesG',0.2,'tStart',0,'onsetOffsetSec',onsetOffsetSec,'verbose',false};
    % point estimate (full data, full perms)
    rPt = groupSeparationTest(N,T,tvec,gargs{:},'nPerms',10000,'seed',1);
    onPt = rPt.onsetTime;
    % bootstrap
    onB = nan(Nboot,1);
    for b=1:Nboot
        iN = randi(nN,nN,1); iT = randi(nT,nT,1);
        rb = groupSeparationTest(N(iN,:),T(iT,:),tvec,gargs{:},'nPerms',nPerms,'seed',b);
        onB(b) = rb.onsetTime;
    end
    detRate = mean(isfinite(onB));
    onBd = onB(isfinite(onB));
    ci = prctile(onBd,[2.5 97.5]);
    results.(meas(m).name) = struct('onsetPoint',onPt,'ci',ci, ...
        'median',median(onBd),'detRate',detRate,'boot',onB, ...
        'nNaive',nN,'nTrained',nT,'tail',tail);
    fprintf('%-9s | onset=%.2fs | 95%% CI [%.2f, %.2f] | bootMed=%.2f | detected %.0f%% of boots | nN=%d nT=%d\n', ...
        meas(m).name, onPt, ci(1), ci(2), median(onBd), 100*detRate, nN, nT);
end

results.tvec = tvec; results.Nboot=Nboot; results.nPerms=nPerms;
save(fullfile(anaDir,'data','onsetCI.mat'),'results');
fprintf('\nSaved onsetCI.mat\n');
