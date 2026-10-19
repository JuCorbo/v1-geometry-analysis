%% computeAxisPositionVsLick.m
%  Figure 6 NoGo-axis-position analysis — cache generator.
%
%  For each TRAINED session, at a fixed post-stimulus snapshot:
%    1. take the main global-z PaCMAP embedding (seed0) at the window nearest 1.5 s
%    2. fit a linear SVM Go(45) vs NoGo -> decision boundary
%    3. signed distance = (X*w + b)/||w||; "corrected" distance sets wrong-side
%       trials to 0 (a NoGo on the Go side -> 0)
%    4. join per-trial lick Response from the public deposit (Trained_behav_slim.csv)
%    5. sort trials by corrected distance, sliding-window bin -> P(lick) vs distance
%
%  RAW only — no Procrustes alignment (the analysis is entirely within-session:
%  per-session SVM, per-session binning, per-session slope — cross-session
%  alignment adds nothing and is dropped).
%
%  Output: cache/DistanceVsLick_<metric>_snap<t>s_bin<b>-<s>.mat  (TallRawTrials, TsumRaw)
clear; clc;

%% ---------------- CONFIG ----------------
seed             = 0;          % main embedding seed
snapshotTime     = 1.5;        % s from stimulus onset
goStim           = 45;
dimsToUse        = [1 2 3];

% Quantified "how Go-like" coordinate to bin lick probability against:
%   'axisProj'     : projection onto the Go->NoGo centroid axis (used in the paper; allows
%                    negatives past the Go centroid; uses stable centroids)
%   'boundaryCorr' : |signed SVM distance| on NoGo side, wrong-side clamped to 0
metric           = 'axisProj';

binSizeTrials    = 20;         % trials per sliding-window point
stepSizeTrials   = 5;          % step between windows
minTrialsPerStim = 10;

svmBoxConstraint = 1;          % SVM kept for the panel-f hyperplane illustration + boundaryCorr
svmStandardize   = false;

%% ---------------- PATHS ----------------
thisDir     = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
figRoot     = fileparts(fileparts(thisDir));                 % .../Figures_final/Figure6
finalRoot   = fileparts(figRoot);                            % .../Figures_final
addpath(fullfile(finalRoot, '_sharedBehavior', 'code'));     % behavDepositTrials

embRoot = fullfile(finalRoot, 'Figure2', 'pacmap', 'cache', 'main_pacmap_bwgz_pm080', ...
                   sprintf('baseline_seed%d', seed), 'BatuhanAngDist');
assert(isfolder(embRoot), 'Embedding root not found:\n  %s', embRoot);

outDir = fullfile(thisDir, '..', 'cache');
if ~exist(outDir, 'dir'), mkdir(outDir); end
outMat = fullfile(outDir, sprintf('DistanceVsLick_%s_snap%.2fs_bin%d-%d.mat', ...
                  metric, snapshotTime, binSizeTrials, stepSizeTrials));

%% ---------------- LICK LOOKUP (deposit) ----------------
% Deposit per-trial Response, Visual block only (SQL-free). The main embeddings
% re-index visual trials 1..N, while the deposit Trial column is continuous across
% the preceding Orientation-Tuning block, so the join is POSITIONAL: embedding
% trial_id t -> t-th deposit visual trial (sorted by Trial). We verify the stim
% sequence matches position-by-position.
B = behavDepositTrials('Trained', true);

%% ---------------- MAIN LOOP ----------------
files = [];
for d = 1:6
    files = [files; dir(fullfile(embRoot, sprintf('D%d', d), '*_PaCMAP_UMATO.mat'))]; %#ok<AGROW>
end
assert(~isempty(files), 'No embedding files under %s', embRoot);

TallRawTrials = table();
TsumRaw       = table();
nUsed = 0; nMissResp = 0; nTot = 0; stimMatch = 0; stimTot = 0;

for fi = 1:numel(files)
    matPath = fullfile(files(fi).folder, files(fi).name);
    S = load(matPath);

    if ~all(isfield(S, {'Y_pacmap','labels','trial_ids','time_window_idx','time_window_center'}))
        warning('Missing fields in %s; skipping.', files(fi).name); continue;
    end
    Y3        = double(S.Y_pacmap(:, dimsToUse));
    labels    = double(S.labels(:));
    trial_ids = double(S.trial_ids(:));
    twi       = double(S.time_window_idx(:));
    twc       = double(S.time_window_center(:));

    experimentID = getExperimentID(S, files(fi).name);
    behavCond    = '';
    if isfield(S,'behavCond'), behavCond = matToChar(S.behavCond); end
    animal = regexp(experimentID, '^[A-Za-z]+\d+', 'match', 'once');

    % Go / NoGo values (expect exactly two stimuli per session export)
    stimVals = unique(labels); stimVals = stimVals(isfinite(stimVals));
    if numel(stimVals) ~= 2 || ~any(stimVals == goStim), continue; end
    goVal   = goStim;
    nogoVal = stimVals(stimVals ~= goStim); nogoVal = nogoVal(1);

    % Snapshot window nearest snapshotTime
    [uW, tC]  = computeWindowCenters(twi, twc);
    [~, k]    = min(abs(tC - snapshotTime));
    snapWIdx  = uW(k); snapT = tC(k);

    idx   = (twi == snapWIdx);
    Xsnap = Y3(idx, :);
    Lsnap = labels(idx);
    Tsnap = trial_ids(idx);
    if ~any(Lsnap == goVal) || ~any(Lsnap == nogoVal), continue; end

    % Per-trial lick Response from deposit (positional join, sorted by Trial)
    d = B(B.Experiment == string(experimentID), :);
    if isempty(d), warning('No deposit trials for %s; skipping.', experimentID); continue; end
    [~, o]  = sort(d.Trial);
    depResp = double(d.Response(o));
    depStim = double(d.Visual_Stim(o));
    Nd = numel(depResp);

    resp = nan(size(Tsnap));
    for i = 1:numel(Tsnap)
        t = round(Tsnap(i));
        if t >= 1 && t <= Nd
            resp(i) = depResp(t);
            stimTot = stimTot + 1;
            stimMatch = stimMatch + (depStim(t) == Lsnap(i));
        end
    end
    nTot = nTot + numel(Tsnap);
    nMissResp = nMissResp + sum(isnan(resp));

    % Fresh per-session linear SVM (boundary illustration + boundaryCorr metric)
    [distAbs, distSigned, ok] = svmAbsDistance(Xsnap, Lsnap, goVal, nogoVal, svmBoxConstraint, svmStandardize);
    if ~ok, continue; end

    % Go->NoGo centroid axis projection (axisProj metric; stable, signed)
    goC   = mean(Xsnap(Lsnap==goVal,  :), 1);
    nogoC = mean(Xsnap(Lsnap==nogoVal,:), 1);
    u     = (nogoC - goC); u = u / norm(u);
    axisProj = (Xsnap - goC) * u';            % 0 at Go centroid; +toward NoGo

    Traw = makeTrialTable(animal, behavCond, experimentID, snapT, ...
                          Tsnap, Lsnap, resp, distAbs, distSigned, axisProj, goVal, nogoVal, metric);
    TallRawTrials = [TallRawTrials; Traw]; %#ok<AGROW>

    Tsum = summarizeBySlidingWindow(Traw, goVal, nogoVal, binSizeTrials, stepSizeTrials, minTrialsPerStim);
    if isempty(Tsum), continue; end
    Tsum.DBname     = repmat("BatuhanAngDist",  height(Tsum), 1);
    Tsum.Animal     = repmat(string(animal),     height(Tsum), 1);
    Tsum.Behav_Cond = repmat(string(behavCond),  height(Tsum), 1);
    Tsum.Experiment = repmat(string(experimentID),height(Tsum), 1);
    Tsum.SnapshotT  = repmat(snapT,              height(Tsum), 1);
    TsumRaw = [TsumRaw; Tsum]; %#ok<AGROW>
    nUsed = nUsed + 1;
end

fprintf('Sessions used: %d / %d files | trials %d | missing Response %d (%.1f%%)\n', ...
    nUsed, numel(files), nTot, nMissResp, 100*nMissResp/max(nTot,1));
fprintf('Stim-sequence consistency (positional join): %.2f%% of %d matched\n', ...
    100*stimMatch/max(stimTot,1), stimTot);

%% ---------------- SAVE ----------------
save(outMat, 'TallRawTrials','TsumRaw','snapshotTime','seed','dimsToUse','metric', ...
    'binSizeTrials','stepSizeTrials','minTrialsPerStim','svmBoxConstraint','svmStandardize','goStim');
fprintf('Saved:\n  %s\n', outMat);

%% ========================= HELPERS =========================
function experimentID = getExperimentID(S, fname)
    experimentID = '';
    if isfield(S,'experimentID') && ~isempty(S.experimentID)
        experimentID = matToChar(S.experimentID);
    end
    if isempty(experimentID)
        [~, base] = fileparts(fname);
        base = regexprep(base, '_PaCMAP_UMATO$', '');
        base = regexprep(base, '^BatuhanAngDist_', '');
        experimentID = base;
    end
end

function out = matToChar(x)
    % Decode IDs stored either as plain char ('C05_121220') or serialized as the
    % text of their ASCII codes ('[67, 48, 53, ...]').
    while iscell(x) && numel(x) == 1, x = x{1}; end
    if isnumeric(x), out = strtrim(char(x(:)')); return; end
    if isstring(x), x = char(x); end
    if ischar(x)
        t = strtrim(x);
        nums = regexp(t, '-?\d+', 'match');
        if ~isempty(nums)
            v = str2double(nums);
            if numel(v) > 1 && all(isfinite(v)) && all((v>=32 & v<=126) | v==9 | v==10 | v==13)
                out = strtrim(char(v(:)')); return;   % ASCII-code list -> text
            end
        end
        out = t; return;
    end
    out = strtrim(char(string(x)));
end

function [uW, tC] = computeWindowCenters(twi, twc)
    uW = sort(unique(twi(isfinite(twi))));
    tC = arrayfun(@(w) mean(twc(twi==w), 'omitnan'), uW);
end

function [distAbs, distSigned, ok] = svmAbsDistance(X, labels, goVal, nogoVal, C, doStd)
    ok = false; distAbs = nan(size(labels)); distSigned = nan(size(labels));
    ybin = nan(size(labels));
    ybin(labels==goVal)   = +1;
    ybin(labels==nogoVal) = -1;
    keep = isfinite(ybin);
    try
        mdl = fitcsvm(X(keep,:), ybin(keep), 'KernelFunction','linear', ...
            'BoxConstraint', C, 'Standardize', doStd, 'ClassNames', [-1 +1]);
        w = mdl.Beta; b = mdl.Bias; den = norm(w);
        if den < eps, return; end
        distSigned = (X*w + b) / den;   % + = Go side
        distAbs    = abs(distSigned);
        ok = true;
    catch
        ok = false;
    end
end

function T = makeTrialTable(animal, behavCond, experimentID, snapT, ...
                            trialIDs, labels, response, distAbs, distSigned, axisProj, goVal, nogoVal, metric)
    n = numel(trialIDs);
    T = table();
    T.Animal      = repmat(string(animal),       n, 1);
    T.Behav_Cond  = repmat(string(behavCond),    n, 1);
    T.Experiment  = repmat(string(experimentID), n, 1);
    T.SnapshotT   = repmat(snapT,                n, 1);
    T.Trial       = trialIDs(:);
    T.Visual_Stim = labels(:);
    T.Response    = response(:);
    T.DistAbs     = distAbs(:);
    T.DistSigned  = distSigned(:);
    T.AxisProj    = axisProj(:);
    % corrected boundary distance: wrong-side trials -> 0
    corr = distAbs(:);
    wrongGo   = (labels(:)==goVal)   & (distSigned(:) < 0);
    wrongNoGo = (labels(:)==nogoVal) & (distSigned(:) > 0);
    corr(wrongGo | wrongNoGo) = 0;
    T.CorrectedDistAbs = corr;
    % Metric = the quantity binned against P(lick) (larger = less Go-like)
    switch metric
        case 'axisProj',     T.Metric = axisProj(:);
        case 'boundaryCorr', T.Metric = corr;
        otherwise, error('Unknown metric %s', metric);
    end
end

function Tsum = summarizeBySlidingWindow(Ttrials, goVal, nogoVal, binSizeTrials, stepSizeTrials, minTrialsPerStim)
    Tsum = table();
    ok = isfinite(Ttrials.Response) & isfinite(Ttrials.Metric);
    T = Ttrials(ok, :);
    if isempty(T), return; end
    for stim = [goVal, nogoVal]
        Ts = T(T.Visual_Stim == stim, :);
        if height(Ts) < minTrialsPerStim, continue; end
        d = Ts.Metric;
        [~, ord] = sort(d, 'ascend');
        n = numel(d);
        if n < binSizeTrials, continue; end
        binIdx = 0;
        for i1 = 1:stepSizeTrials:(n - binSizeTrials + 1)
            winOrd = ord(i1 : i1 + binSizeTrials - 1);
            binIdx = binIdx + 1;
            dd = Ts.Metric(winOrd);
            rr = table();
            rr.Stim       = stim;
            rr.Bin        = binIdx;
            rr.nTrials    = numel(winOrd);
            rr.lickProb   = mean(Ts.Response(winOrd) == 1, 'omitnan');
            rr.distMedian = median(dd, 'omitnan');
            rr.distMin    = min(dd, [], 'omitnan');
            rr.distMax    = max(dd, [], 'omitnan');
            Tsum = [Tsum; rr]; %#ok<AGROW>
        end
    end
end
