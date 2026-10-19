%% pull_pooledSpectralPR.m  (companion to pull_SVMPCATimeCourse_updated.m)
% Pooled-classes ("overall") participation ratio: PCA on the Go and NoGo trials
% POOLED and centred on the pooled mean, so the spectrum contains both the
% within-class variability and the between-class separation.
% Same data path and subsampling as pull_SVMPCATimeCourse_updated.m:
% getWindowed('baseline') (sigma floor 0.12), 250-cell ceiling, 100 trials per
% class, 5 rounds with the SAME seeds (rng(1000*rep+10*day+sess) for cells,
% rng(100*rep+10*day+sess) for trials) -> identical subsamples to dimMetrics.
% Saves cache/pooledSpectralPR.mat: PRpool{dbi,day}{sess} = [nRounds x nW].
clear; clc;
codeDir    = fileparts(mfilename('fullpath'));
aDir       = fileparts(codeDir);                       % dimensionality
figDir     = fileparts(fileparts(aDir));               % Figures_final
addpath(fullfile(figDir,'_sharedWindowedData','code'));
sharedData = fullfile(figDir,'_sharedWindowedData','data');

DBlist   = {'BatuhanADNaive','BatuhanAngDist'};        % naive(1), trained(2)
groupDir = {'naive','trained'};
BehavCond= {'D1','D2','D3','D4','D5','D6'};
GoOri = 45;
nRounds = 5; nCells = 250; nSamplesMax = 100;

PRpool = cell(2, numel(BehavCond)); sessNames = cell(2, numel(BehavCond));
tvec = [];
for dbi = 1:2
    for exps = 1:numel(BehavCond)
        files = dir(fullfile(sharedData, groupDir{dbi}, sprintf('D%d',exps), '*_win.mat'));
        PRpool{dbi,exps} = {}; sessNames{dbi,exps} = {};
        for sess = 1:numel(files)
            S = load(fullfile(files(sess).folder, files(sess).name));
            [X, info] = getWindowed(S, 'baseline');     % [nTr x nC x nW]
            if isempty(tvec), tvec = info.tvec(:)'; end
            vs = double(info.Visual_Stim(:));
            go = find(vs == GoOri); nogo = find(vs ~= GoOri & isfinite(vs));
            if numel(go) < 2 || numel(nogo) < 2, continue; end
            [~, nC, nW] = size(X); nUseCells = min(nC, nCells);
            PRrep = NaN(nRounds, nW);
            for repet = 1:nRounds
                rng(1000*repet + 10*exps + sess);
                cellSub = randperm(nC, nUseCells);
                nSamples = min([nSamplesMax, numel(go), numel(nogo)]);
                rng(100*repet + 10*exps + sess);
                idx1 = randsample(go, nSamples); idx2 = randsample(nogo, nSamples);
                for w = 1:nW
                    Xw = X(:, cellSub, w);
                    Xp = [Xw(idx1,:); Xw(idx2,:)];               % pooled classes
                    Xp = Xp - mean(Xp, 1, 'omitnan');            % pooled mean centring
                    lam = svd(double(Xp), 'econ').^2;            % eigenvalues of X'X (prop. to covariance spectrum)
                    lam(lam < 0) = 0; s1 = sum(lam); s2 = sum(lam.^2);
                    if s2 > 0, PRrep(repet, w) = (s1^2)/s2; end
                end
            end
            PRpool{dbi,exps}{end+1} = PRrep;
            sessNames{dbi,exps}{end+1} = files(sess).name;
        end
        fprintf('  %-7s D%d: %d sessions\n', groupDir{dbi}, exps, numel(PRpool{dbi,exps}));
    end
end

outMat = fullfile(aDir,'cache','pooledSpectralPR.mat');
save(outMat, 'PRpool','sessNames','tvec','BehavCond','DBlist','-v7.3');
fprintf('saved %s\n', outMat);
