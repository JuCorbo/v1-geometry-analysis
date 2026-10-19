function pull_forDR_fromCache(scheme, outRoot)
%PULL_FORDR_FROMCACHE  Build the PaCMAP/UMATO input files from the shared cache.
% For every shared-cache session it calls
% getWindowed(scheme) and stacks the 26 windows into the [nTrials*nWin x nCells]
% row layout the Python embedder (s02_make_pacmap_umato.py) and the MATLAB renderers
% expect, writing <DB>_<session>_TimeDR_forDR.mat under outRoot/<DB>/D<d>/.
%
%   scheme   getWindowed scheme for X ('global' main, 'baseline' supp/multiDR)
%   outRoot  destination root (per-DB / per-day subfolders created here)
%
% Row order: for w=1..nWin, append the [nTrials x nCells]
% block; labels=Visual_Stim, trial_ids=1..nTrials, time_window_idx=w,
% time_window_center=tvec(w) (window last-frame time).

codeDir    = fileparts(mfilename('fullpath'));
figDir     = fileparts(fileparts(fileparts(codeDir)));      % Figures_final
sharedCode = fullfile(figDir,'_sharedWindowedData','code');
sharedData = fullfile(figDir,'_sharedWindowedData','data');
addpath(sharedCode);

files = dir(fullfile(sharedData,'*','D*','*_win.mat'));
fprintf('pull_forDR_fromCache(scheme=%s): %d sessions -> %s\n', scheme, numel(files), outRoot);

for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name));
    [X,  info] = getWindowed(S, scheme);        % [nTr x nC x nW]
    [nTr, nC, nW] = size(X);

    X_all   = zeros(nTr*nW, nC);
    labels             = zeros(nTr*nW, 1);
    trial_ids          = zeros(nTr*nW, 1);
    time_window_idx    = zeros(nTr*nW, 1);
    time_window_center = zeros(nTr*nW, 1);
    vs = double(info.Visual_Stim(:));
    tv = double(info.tvec(:));
    for w = 1:nW
        r = (w-1)*nTr + (1:nTr);
        X_all(r,:) = X(:,:,w);
        labels(r)             = vs;
        trial_ids(r)          = (1:nTr)';
        time_window_idx(r)    = w;
        time_window_center(r) = tv(w);
    end

    % --- metadata (passthrough for renderers) ---
    DBname       = char(info.DB);
    experimentID = char(info.session);
    tok = regexp(files(k).folder, 'D(\d)$', 'tokens', 'once');
    dnum = str2double(tok{1});
    behavCond    = sprintf('D%d', dnum);
    Vstim_all    = vs;

    % --- save as the Python/renderers expect ---
    base    = erase(files(k).name, '_win.mat');            % <DB>_<session>
    dayDir  = fullfile(outRoot, DBname, sprintf('D%d', dnum));
    if ~exist(dayDir,'dir'), mkdir(dayDir); end
    outPath = fullfile(dayDir, [base '_TimeDR_forDR.mat']);

    X = X_all;   %#ok<NASGU>  save under the name the embedder reads
    save(outPath, 'X', 'labels', 'trial_ids', ...
        'time_window_idx', 'time_window_center', 'Vstim_all', ...
        'DBname', 'experimentID', 'behavCond', '-v7.3');
end
fprintf('  done (%s).\n', scheme);
end
