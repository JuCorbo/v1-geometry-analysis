%% locoCapacityCompute.m  — PHASE 2 of 2 (capacity; NO DB)
% Loads each locoCapacityPrep.m output and runs check_capacity for the T1/T3/EQ
% partitions across the rolling-window schedule, in a clean process (no DB, no
% signal toolbox, no byte-stream decode) -- which is the ONLY context in which
% check_capacity does not segfault (see locoCapacityPrep.m header).
%
% Capacity logic is VERBATIM from the lick/pupil pipeline.
% Output: D:/embedding/TraceSpikes/D#/locoCapacity/<sess>_<DB>_post_<part>_capacity_timecourse.mat
% Resumable: skips partitions whose output .mat already exists.
%
% Overrides: RUN_PREP_ONLY (name substring to select one prep file),
%            RUN_CENTERS_SUBSET, RUN_RESUME, RUN_SAVEROOT (keep SHORT paths).

clearvars -except RUN_PREP_ONLY RUN_CENTERS_SUBSET RUN_RESUME RUN_SAVEROOT; clc; close all;

if exist('check_capacity','file') == 0
    capLibRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))),'toolbox_dnn-object-manifolds');
    if exist(capLibRoot,'dir'), addpath(genpath(capLibRoot)); end
    assert(exist('check_capacity','file')>0, 'check_capacity not found (%s)', capLibRoot);
end
if exist('cplexoptimset','file') == 0
    cplexDir = 'C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64';
    if exist(cplexDir,'dir'), addpath(cplexDir); end
    assert(exist('cplexoptimset','file')>0, 'cplexoptimset (CPLEX) not found (%s)', cplexDir);
end

%% ---------------- SETTINGS ----------------
DB = 'BatuhanAngDist';
prepRoot = 'D:/embedding/TraceSpikes/locoPrep';
capacitySuffix = '_capacity_balanced_s1.mat';
minTrialsPerCondPart = 8;

saveRoot = 'D:/embedding/TraceSpikes/';            % SHORT paths (do not use long scratch)
if exist('RUN_SAVEROOT','var') && ~isempty(RUN_SAVEROOT), saveRoot = RUN_SAVEROOT; end
tmpRoot  = fullfile(saveRoot, 'locoCapTmp');
if ~exist(tmpRoot,'dir'), mkdir(tmpRoot); end

if ~exist('RUN_PREP_ONLY','var'),     RUN_PREP_ONLY = ''; end
if ~exist('RUN_CENTERS_SUBSET','var'),RUN_CENTERS_SUBSET = []; end
if ~exist('RUN_RESUME','var'),        RUN_RESUME = true; end

prepFiles = dir(fullfile(prepRoot, sprintf('*_%s_prep.mat', DB)));
if ~isempty(RUN_PREP_ONLY)
    prepFiles = prepFiles(contains({prepFiles.name}, RUN_PREP_ONLY));
end
fprintf('Found %d prep files.\n', numel(prepFiles));

pnames = {'T1','T3','EQ'};
for f = 1:numel(prepFiles)
    P = load(fullfile(prepFiles(f).folder, prepFiles(f).name));
    sess = P.sess; day = P.day; goVal = P.goVal;
    centers = P.centers; rollingWindows = P.rollingWindows; tvec = P.tvec;
    if ~isempty(RUN_CENTERS_SUBSET)
        centers = centers(RUN_CENTERS_SUBSET); rollingWindows = rollingWindows(RUN_CENTERS_SUBSET); tvec = tvec(RUN_CENTERS_SUBSET);
    end
    fprintf('\n=== COMPUTE %s | %s ===\n', day, sess);

    saveDir = fullfile(saveRoot, ['D' day(end)], 'locoCapacity');
    if ~exist(saveDir,'dir'), mkdir(saveDir); end
    partSets = struct('T1',P.T1,'T3',P.T3,'EQ',P.EQ);

    for pi = 1:numel(pnames)
        pn = pnames{pi}; sel = partSets.(pn);
        nGo = sum(P.isGo(sel)); nNo = sum(~P.isGo(sel));
        if nGo < minTrialsPerCondPart || nNo < minTrialsPerCondPart
            warning('%s | %s: too few (Go=%d NoGo=%d). Skipping.', sess, pn, nGo, nNo); continue;
        end
        outFile = fullfile(saveDir, sprintf('%s_%s_loco_post_%s_capacity_timecourse.mat', sess, DB, pn));
        if RUN_RESUME && isfile(outFile)
            fprintf('  resume: skip existing %s\n', pn); continue;
        end
        Tsel = P.T(sel,:,:); Vsel = P.Vstim(sel);
        cap = NaN(numel(centers),1);
        for k = 1:numel(centers)
            prefix = fullfile(tmpRoot, sprintf('%s_%s_post_%s_t%03d', sess, DB, pn, centers(k)));
            cap(k) = capacity_from_window(Tsel, Vsel, rollingWindows{k}, prefix, capacitySuffix, goVal);
        end
        capacity_rolling = cap; %#ok<NASGU>
        meta = P.meta; meta.partition = pn; meta.nGo = nGo; meta.nNoGo = nNo; %#ok<NASGU>
        save(outFile,'capacity_rolling','tvec','meta','-v7.3');
        fprintf('  saved %s (Go=%d NoGo=%d)\n', pn, nGo, nNo);
    end
end
fprintf('\nCompute done.\n');

%% ================= LOCAL FUNCTIONS =================
function cap = capacity_from_window(T, Vstim, winIdx, name_prefix, capacity_suffix, goVal)
% VERBATIM capacity logic (lick/pupil pipeline).
    cap = NaN;
    if isempty(winIdx), return; end
    X = squeeze(nansum(T(:,:,winIdx),3));
    if isvector(X), X = reshape(X, numel(Vstim), []); end
    isGo = (Vstim==goVal); isNoGo = ~isGo;
    nG = sum(isGo); nN = sum(isNoGo);
    if nG==0 || nN==0, return; end
    if nG >= nN
        TF = cat(3, X(isGo,:)', [X(isNoGo,:); NaN(nG-nN, size(X,2))]');
    else
        TF = cat(3, [X(isGo,:); NaN(nN-nG, size(X,2))]', X(isNoGo,:)');
    end
    tuning_function = TF; %#ok<NASGU>
    data_titles = {name_prefix}; %#ok<NASGU>
    session_ids = 1;
    tuneFile = [name_prefix '_tuning.mat']; capFile = [name_prefix capacity_suffix];
    if isfile(tuneFile), delete(tuneFile); end
    if isfile(capFile),  delete(capFile);  end
    save(tuneFile, 'tuning_function', 'data_titles', '-v7.3');
    try
        check_capacity(name_prefix, session_ids);
    catch ME
        warning('check_capacity failed at %s: %s', name_prefix, ME.message); return;
    end
    if ~isfile(capFile)
        d = dir([fileparts(name_prefix) filesep '*' capacity_suffix]);
        if isempty(d), return; else, capFile = fullfile(d(1).folder, d(1).name); end
    end
    file = matfile(capFile);
    cap = 2 ./ file.capacity_results;
end
