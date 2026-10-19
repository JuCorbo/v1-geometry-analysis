%% pull_rollingCapacity.m
% Manifold-capacity ROLLING-window timecourse per session, computed from the
% shared windowed intermediate.
%
% Data source: getWindowed(S,'sub') = Xr - 3*b (per-trial baseline-subtracted
% 3-frame window sums; summing a baseline-subtracted trace over a 3-frame
% window equals (raw 3-frame sum) - 3*(baseline mean) = Xr - 3*b). The rolling
% centres are 9:3:84, exactly the 26 windows of the shared cache.
%
% NOTE ON RUNNING THIS: the deliverable figures read the 324 caches shipped in
% ../cache/capacityTimecourse; this script can REGENERATE them from
% getWindowed('sub') plus the bundled Chung/Sompolinsky toolbox
% (_capacity_source/toolbox_dnn-object-manifolds). check_capacity reads only
% the *_tuning.mat this script writes; the one external dependency is a SOLVER:
% the toolbox's separability test (check_linear_seperability_svm_cplexqp)
% requires IBM CPLEX (cplexqp/cplexoptimset). Without CPLEX it errors at the
% separability step. By default it writes to a '_regen' scratch folder so a
% stray run cannot overwrite the shipped caches; set CAP_SAVEROOT to
% regenerate for real.
clear; clc;

%% ---------------------- PATHS (shared cache) ----------------------
codeDir    = fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir = pwd; end
anaDir     = fileparts(codeDir);                        % manifoldCapacity
figDir     = fileparts(fileparts(anaDir));              % Figures_final
sharedCode = fullfile(figDir,'_sharedWindowedData','code');
sharedData = fullfile(figDir,'_sharedWindowedData','data');
addpath(sharedCode);
addpath(genpath(fullfile(anaDir,'_capacity_source','toolbox_dnn-object-manifolds')));  % check_capacity + Chung/Sompolinsky replica-MFT toolbox (bundled)
addpath(codeDir);

%% ---------------------- USER PARAMS ----------------------
DBlist   = {'BatuhanADNaive','BatuhanAngDist'};   % naive(1), trained(2)
groupDir = {'naive','trained'};
BehavCond= {'D1','D2','D3','D4','D5','D6'};

% Output root: default to a scratch '_regen' folder so we never clobber the
% retained 324 caches. Override with CAP_SAVEROOT to regenerate for real.
saveRoot = getenv('CAP_SAVEROOT');
if isempty(saveRoot), saveRoot = fullfile(anaDir,'capacityTimecourse_regen'); end

capacity_suffix = '_capacity_balanced_s1.mat';   % exact filename check_capacity writes
GoOri = 45;

%% ---------------------- PER GROUP / DAY / SESSION LOOPS ----------------------
for dbi = 1:2
    DB = DBlist{dbi};
    for exps = 1:numel(BehavCond)
        tic;
        files = dir(fullfile(sharedData, groupDir{dbi}, sprintf('D%d',exps), '*_win.mat'));
        if isempty(files), warning('No sessions for %s %s', groupDir{dbi}, BehavCond{exps}); continue; end

        saveFolderDay = fullfile(saveRoot, ['D' num2str(exps)], 'timecourse');
        if ~exist(saveFolderDay, 'dir'), mkdir(saveFolderDay); end

        for sess = 1:numel(files)
            S = load(fullfile(files(sess).folder, files(sess).name));
            [X, info] = getWindowed(S, 'sub');        % [nTrials x nCells x nWin] = Xr - 3*b
            sessName  = char(info.session);
            Vstim     = double(info.Visual_Stim(:));
            tvec      = info.tvec(:)';
            nW        = size(X,3);

            fprintf('\n=== %s D%d  Session %d/%d (%s) ===\n', ...
                groupDir{dbi}, exps, sess, numel(files), sessName);

            % Need both classes present
            if ~any(Vstim==GoOri) || ~any(Vstim~=GoOri & isfinite(Vstim))
                warning('Skipping %s: missing Go or NoGo trials.', sessName); continue;
            end

            % ----------- Compute rolling capacity, one value per window -----------
            capacity_rolling = NaN(nW,1);
            for k = 1:nW
                Xw = X(:,:,k);                        % [nTrials x nCells], already the window sum
                name_roll = fullfile(saveFolderDay, sprintf('%s_%s_roll_w%03d', sessName, DB, k));
                capacity_rolling(k) = capacity_from_windowX(Xw, Vstim, GoOri, name_roll, capacity_suffix);
            end

            % ----------- Save summary (schema matches the retained caches) -----------
            meta = struct('DB',DB,'session',sessName,'BehavCond',BehavCond{exps}, ...
                          'dt',info.dt,'subtractBaseline',true, ...
                          'winFrames',info.winFrames,'stepFrames',3, ...
                          'source','getWindowed(sub) shared cache (SQL-free)');

            outFile = fullfile(saveFolderDay,[sessName '_' DB '_capacity_timecourse.mat']);
            save(outFile,'capacity_rolling','tvec','meta','-v7.3');
            fprintf('Saved %s\n', outFile);
        end
        fprintf('Finished %s %s in %.1f s\n', groupDir{dbi}, BehavCond{exps}, toc);
    end
end
disp('All done.');

%% ---------------------- LOCAL FUNCTIONS ----------------------
function cap = capacity_from_windowX(X, Vstim, GoOri, name_prefix, capacity_suffix)
    % X is already the [nTrials x nCells] baseline-subtracted window sum.
    % Balance Go/NoGo into tuning_function [nCells x nTrials x 2] and run capacity.
    cap = NaN;
    isGo = (Vstim==GoOri); isNoGo = ~isGo & isfinite(Vstim);
    nG = sum(isGo); nN = sum(isNoGo);
    if nG==0 || nN==0, return; end

    if nG >= nN
        TF = cat(3, X(isGo,:)', [X(isNoGo,:); NaN(nG-nN, size(X,2))]');
    else
        TF = cat(3, [X(isGo,:); NaN(nN-nG, size(X,2))]', X(isNoGo,:)');
    end
    tuning_function = TF; %#ok<NASGU>
    data_titles     = {name_prefix}; %#ok<NASGU>
    session_ids     = 1;

    save([name_prefix '_tuning.mat'], 'tuning_function', 'data_titles', '-v7.3');
    try
        check_capacity(name_prefix, session_ids);
    catch ME
        warning('check_capacity failed at %s: %s', name_prefix, ME.message);
        return;
    end

    filePath = [name_prefix capacity_suffix];
    if ~isfile(filePath)
        d = dir([fileparts(name_prefix) filesep '*' capacity_suffix]);
        if isempty(d), return; else, filePath = fullfile(d(1).folder, d(1).name); end
    end
    try
        mf = matfile(filePath);
        cap = 2 ./ mf.capacity_results;
        if numel(cap) > 1, cap = cap(1); end
    catch ME
        warning('Failed reading capacity_results in %s: %s', filePath, ME.message);
        cap = NaN;
    end
end
