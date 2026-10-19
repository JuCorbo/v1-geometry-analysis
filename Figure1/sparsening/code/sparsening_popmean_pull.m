%% sparsening_popmean_pull.m  build the slim pop-mean intermediate from the shared cache
% Reads the shared windowed intermediate via getWindowed('baseline' = baseline
% z-score, sigma floored at 0.12) and produces the slim cache the make script
% consumes (per-session pop-mean timecourses + two example sessions).
% Population mean = mean across neurons of the trial-averaged, baseline-z, window-integrated
% signal, per Go/NoGo, per session. Plot axis: BOTH tvec (last frame) and tvec_centre.
codeDir = fileparts(mfilename('fullpath'));
aDir    = fileparts(codeDir);
figDir  = fileparts(fileparts(aDir));                 % Figures_final
sharedCode = fullfile(figDir,'_sharedWindowedData','code');
sharedData = fullfile(figDir,'_sharedWindowedData','data');
addpath(sharedCode);
cacheDir = fullfile(aDir,'cache'); if ~exist(cacheDir,'dir'), mkdir(cacheDir); end

exampleNaiveFile = fullfile(sharedData,'naive','D3','BatuhanADNaive_056_20623_0_win.mat');
exampleTrainFile = fullfile(sharedData,'trained','D3','BatuhanAngDist_F07_42221__win.mat');

groups = {'naive','trained'};
allGo  = {[],[]};   % {naive, trained} : [nSess x nWin]
allNo  = {[],[]};
tvec = []; tvec_centre = [];

for gi = 1:2
    files = dir(fullfile(sharedData,groups{gi},'D*','*_win.mat'));
    for k = 1:numel(files)
        S = load(fullfile(files(k).folder,files(k).name));
        [X, info] = getWindowed(S,'baseline');        % [nTrials x nCells x nWin]
        if isempty(tvec), tvec = info.tvec(:)'; tvec_centre = info.tvec_centre(:)'; end
        vs  = info.Visual_Stim(:);
        goM = (vs==45); noM = ~goM;                   % 2 stimuli per session
        tcGo = squeeze(mean(mean(X(goM,:,:),1,'omitnan'),2,'omitnan'))';   % 1 x nWin (neurons then trials)
        tcNo = squeeze(mean(mean(X(noM,:,:),1,'omitnan'),2,'omitnan'))';
        allGo{gi}(end+1,:) = tcGo;    %#ok<SAGROW>
        allNo{gi}(end+1,:) = tcNo;    %#ok<SAGROW>
    end
    fprintf('  %-7s: %d sessions\n', groups{gi}, size(allGo{gi},1));
end
allNaiveGo = allGo{1}; allNaiveNo = allNo{1};
allTrainGo = allGo{2}; allTrainNo = allNo{2};

% ---- example sessions: per-neuron trial-mean rasters + pop-mean TC +/- SEM(over neurons) ----
exampleNaive = buildExample(exampleNaiveFile,'056_20623_0');
exampleTrain = buildExample(exampleTrainFile,'F07_42221');

save(fullfile(cacheDir,'sparsening_popmean_pulled.mat'), ...
    'allNaiveGo','allNaiveNo','allTrainGo','allTrainNo', ...
    'exampleNaive','exampleTrain','tvec','tvec_centre','-v7.3');
fprintf('saved slim pop-mean intermediate -> %s\n', fullfile(cacheDir,'sparsening_popmean_pulled.mat'));

function ex = buildExample(fpath,name)
    S = load(fpath); [X, info] = getWindowed(S,'baseline');
    vs = info.Visual_Stim(:); goM = (vs==45); noM = ~goM;
    matGo = squeeze(mean(X(goM,:,:),1,'omitnan'));    % nCells x nWin (per-neuron trial mean)
    matNo = squeeze(mean(X(noM,:,:),1,'omitnan'));
    tc = info.tvec(:)'; early = tc>=0 & tc<1;         % sort neurons by early (0-1 s) Go activity, desc
    [~,ord] = sort(mean(matGo(:,early),2,'omitnan'),'descend');
    matGo = matGo(ord,:); matNo = matNo(ord,:);
    semN  = @(M) std(M,0,1,'omitnan')./sqrt(sum(~isnan(M),1));
    ex = struct('name',name, 'matGo',matGo,'matNo',matNo, ...
        'tcGo',mean(matGo,1,'omitnan'),'semGo',semN(matGo), ...
        'tcNo',mean(matNo,1,'omitnan'),'semNo',semN(matNo));
end
