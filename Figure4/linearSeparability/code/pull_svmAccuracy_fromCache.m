%% pull_svmAccuracy_fromCache.m
% Full-data linear-SVM accuracy timecourse (the array make_svmFullTimecourse.m
% consumes), computed from the shared windowed cache via getWindowed('baseline')
% (baseline z-score, sigma floor 0.12).
%
% Decoding: Go(45) vs the day's NoGo orientation, up to 250 cells subsampled per
% session, up to 100 trials/class balanced, nRounds=5, linear SVM, 5-fold
% cross-validation, accuracy in percent per fold. Same per-session/round rng
% seeds as the Figure-3 dimensionality pull, so the cell/trial subsamples are
% identical between the two pulls.
%
% Saves accuracy [nDB x nDays x nSess x nRounds x nWin x nCV] + tvec to cache/.
clear; clc;
codeDir    = fileparts(mfilename('fullpath'));
anaDir     = fileparts(codeDir);                        % linearSeparability
figDir     = fileparts(fileparts(anaDir));              % Figures_final
sharedCode = fullfile(figDir,'_sharedWindowedData','code');
sharedData = fullfile(figDir,'_sharedWindowedData','data');
addpath(sharedCode);

DBlist   = {'BatuhanADNaive','BatuhanAngDist'};   % naive(1), trained(2)
groupDir = {'naive','trained'};
BehavCond= {'D1','D2','D3','D4','D5','D6'};
GoOri = 45;
nRounds = 5; nCells = 250; nSamplesMax = 100; nCV = 5;

% Discover the window count / tvec from the first available session.
tvec = []; nWin = [];
probe = dir(fullfile(sharedData, groupDir{1}, 'D1', '*_win.mat'));
S0 = load(fullfile(probe(1).folder, probe(1).name));
[~, i0] = getWindowed(S0, 'baseline'); tvec = i0.tvec(:)'; nWin = numel(tvec);

maxSess  = 0;
for dbi = 1:2, for exps = 1:6
    maxSess = max(maxSess, numel(dir(fullfile(sharedData, groupDir{dbi}, sprintf('D%d',exps), '*_win.mat'))));
end, end

accuracy = NaN(2, 6, maxSess, nRounds, nWin, nCV);

for dbi = 1:2
    for exps = 1:numel(BehavCond)
        files = dir(fullfile(sharedData, groupDir{dbi}, sprintf('D%d',exps), '*_win.mat'));
        for sess = 1:numel(files)
            S = load(fullfile(files(sess).folder, files(sess).name));
            [X, info] = getWindowed(S, 'baseline');    % [nTr x nC x nW], baseline-z @ 0.12
            vs = double(info.Visual_Stim(:));
            go   = find(vs == GoOri);
            nogo = find(vs ~= GoOri & isfinite(vs));
            if numel(go) < 2 || numel(nogo) < 2, continue; end
            [~, nC, ~] = size(X); nUseCells = min(nC, nCells);
            for repet = 1:nRounds
                rng(1000*repet + 10*exps + sess);
                cellSub = randperm(nC, nUseCells);
                nSamples = min([nSamplesMax, numel(go), numel(nogo)]);
                rng(100*repet + 10*exps + sess);
                idx1 = randsample(go, nSamples); idx2 = randsample(nogo, nSamples);
                targetVec = [ones(nSamples,1); zeros(nSamples,1)];
                for w = 1:nWin
                    Xw = X(:, cellSub, w);
                    data_win = [Xw(idx1,:); Xw(idx2,:)];
                    if all(~isfinite(data_win(:))) || all(data_win(:) == 0), continue; end
                    CVSVMModel = fitcsvm(data_win, targetVec, ...
                        'KernelFunction','linear', 'CrossVal','on', 'KFold', nCV);
                    for cv = 1:nCV
                        testIdx = CVSVMModel.Partition.test(cv);
                        predLab = predict(CVSVMModel.Trained{cv}, data_win(testIdx,:));
                        accuracy(dbi, exps, sess, repet, w, cv) = mean(predLab == targetVec(testIdx)) * 100;
                    end
                end
            end
        end
        fprintf('  %-7s D%d: %d sessions done\n', groupDir{dbi}, exps, numel(files));
    end
end

outMat = fullfile(anaDir,'cache','svm_accuracy_baseline012.mat');
if ~exist(fullfile(anaDir,'cache'),'dir'), mkdir(fullfile(anaDir,'cache')); end
save(outMat, 'accuracy','tvec','BehavCond','DBlist','-v7.3');
fprintf('saved %s\n', outMat);
