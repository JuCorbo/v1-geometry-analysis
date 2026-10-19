%% pull_SVMPCATimeCourse_updated.m
% Dimensionality pull: computes dimMetrics (per-class PCA -> D_PR, k80, tail
% slope, per session/day/round/window/class) from the shared windowed cache via
% getWindowed('baseline') (baseline z-score, sigma floor 0.12).
% Subsampling: nCells 250 ceiling, nSamples 100/class, nRounds 5, seeded per
% session/round (same seeds as the Fig 4 SVM pull). Saves dimMetrics to cache/.
clear; clc;
codeDir    = fileparts(mfilename('fullpath'));
aDir       = fileparts(codeDir);                       % dimensionality
figDir     = fileparts(fileparts(aDir));               % Figures_final
sharedCode = fullfile(figDir,'_sharedWindowedData','code');
sharedData = fullfile(figDir,'_sharedWindowedData','data');
addpath(sharedCode);

DBlist   = {'BatuhanADNaive','BatuhanAngDist'};        % naive(1), trained(2)
groupDir = {'naive','trained'};
BehavCond= {'D1','D2','D3','D4','D5','D6'};
GoOri = 45;
nRounds = 5; nCells = 250; nSamplesMax = 100;
maxPCSave = 80; tailFitStart = 10; tailFitEnd = 40;

dimMetrics = struct('go',{},'nogo',{},'nUseCells',{},'nSamples',{},'N_D_ratio',{},'tvec',{});
tvec = [];

for dbi = 1:2
    for exps = 1:numel(BehavCond)
        files = dir(fullfile(sharedData, groupDir{dbi}, sprintf('D%d',exps), '*_win.mat'));
        for sess = 1:numel(files)
            S = load(fullfile(files(sess).folder, files(sess).name));
            [X, info] = getWindowed(S, 'baseline');     % [nTr x nC x nW]
            if isempty(tvec), tvec = info.tvec(:)'; end
            vs = double(info.Visual_Stim(:));
            go = find(vs == GoOri); nogo = find(vs ~= GoOri & isfinite(vs));
            if numel(go) < 2 || numel(nogo) < 2, continue; end
            [nTr, nC, nW] = size(X); nUseCells = min(nC, nCells);
            for repet = 1:nRounds
                rng(1000*repet + 10*exps + sess);
                cellSub = randperm(nC, nUseCells);
                nSamples = min([nSamplesMax, numel(go), numel(nogo)]);
                rng(100*repet + 10*exps + sess);
                idx1 = randsample(go, nSamples); idx2 = randsample(nogo, nSamples);
                dim_go_tc   = repmat(empty_dim_struct(maxPCSave,nSamples,nUseCells), nW, 1);
                dim_nogo_tc = repmat(empty_dim_struct(maxPCSave,nSamples,nUseCells), nW, 1);
                for w = 1:nW
                    Xw = X(:, cellSub, w);
                    Xg = Xw(idx1,:); Xn = Xw(idx2,:);
                    Xg = Xg - mean(Xg,1,'omitnan'); Xn = Xn - mean(Xn,1,'omitnan');
                    dim_go_tc(w)   = compute_dim_pca(Xg, maxPCSave, tailFitStart, tailFitEnd);
                    dim_nogo_tc(w) = compute_dim_pca(Xn, maxPCSave, tailFitStart, tailFitEnd);
                end
                dimMetrics(dbi,exps,sess,repet).go = dim_go_tc;
                dimMetrics(dbi,exps,sess,repet).nogo = dim_nogo_tc;
                dimMetrics(dbi,exps,sess,repet).nUseCells = nUseCells;
                dimMetrics(dbi,exps,sess,repet).nSamples  = nSamples;
                dimMetrics(dbi,exps,sess,repet).N_D_ratio = nSamples / nUseCells;
                dimMetrics(dbi,exps,sess,repet).tvec = tvec;
            end
        end
        fprintf('  %-7s D%d: %d sessions\n', groupDir{dbi}, exps, numel(files));
    end
end

outMat = fullfile(aDir,'cache','dimMetrics_baseline012.mat');
if ~exist(fullfile(aDir,'cache'),'dir'), mkdir(fullfile(aDir,'cache')); end
save(outMat, 'dimMetrics','tvec','BehavCond','DBlist','-v7.3');
fprintf('saved %s\n', outMat);

%% ---- PCA helpers (from SVMPCATimeCourse_updated.m) ----
function out = compute_dim_pca(X, maxPCSave, tailFitStart, tailFitEnd)
X = double(X); N = size(X,1); D = size(X,2);
if N < 3 || D < 2, out = empty_dim_struct(maxPCSave, N, D); return; end
[~, ~, latent, ~, explained] = pca(X, 'Centered', false);
P = numel(latent); Puse = min(P, maxPCSave);
latPad = NaN(1,maxPCSave); expPad = NaN(1,maxPCSave);
latPad(1:Puse) = latent(1:Puse); expPad(1:Puse) = explained(1:Puse);
cumEV = cumsum(explained);
k80 = find(cumEV >= 80, 1, 'first'); if isempty(k80), k80 = NaN; end
lam = latent(:); lam(lam < 0) = 0; s1 = sum(lam); s2 = sum(lam.^2);
D_PR = NaN; if s2 > 0, D_PR = (s1^2)/s2; end
tailSlope = NaN; a = max(2,tailFitStart); b = min(P,tailFitEnd);
if b > a, idx = (a:b)'; yy = log(max(lam(a:b),1e-12)); pp = polyfit(idx,yy,1); tailSlope = pp(1); end
out.latent=latPad; out.explained=expPad; out.k80=k80; out.D_PR=D_PR;
out.tailSlope=tailSlope; out.Pmax=P; out.N=N; out.D=D; out.N_D_ratio=N/D;
end

function out = empty_dim_struct(maxPCSave, N, D)
out.latent=NaN(1,maxPCSave); out.explained=NaN(1,maxPCSave);
out.k80=NaN; out.D_PR=NaN; out.tailSlope=NaN; out.Pmax=NaN;
out.N=N; out.D=D; out.N_D_ratio=N/max(D,1);
end
