%% computeNativeAxisVsLick.m
%  Figure 6 NATIVE-SPACE supplement (identical pipeline to the PaCMAP axis-
%  position version, but the Go->NoGo axis is computed in the full native
%  neural space rather than in a 3-D embedding).
%
%  For each TRAINED session, at the window nearest 1.5 s:
%    1. native population activity X = getWindowed(session,'global')  [nTrials x nCells]
%    2. Go->NoGo centroid axis u = (nogoC - goC)/||.||   (in nCells-dim space)
%    3. axisProj = (x - goC)*u   per trial (Near=Go-like -> Far)
%    4. lick Response comes straight from the cache (info.Response), no join needed
%    5. sort NoGo trials by axisProj, sliding-window bin -> P(lick) vs axis position
%
%  Output: cache/NativeAxisVsLick_global_snap1.50s_bin20-5.mat  (TallRawTrials, TsumRaw)
clear; clc;

%% ---------------- CONFIG ----------------
scheme           = 'global';   % same normalisation the PaCMAP main was built from
snapshotTime     = 1.5;
goStim           = 45;
binSizeTrials    = 20;
stepSizeTrials   = 5;
minTrialsPerStim = 10;

thisDir   = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
figRoot   = fileparts(fileparts(thisDir));                 % .../Figures_final/Figure6
finalRoot = fileparts(figRoot);                            % .../Figures_final
addpath(fullfile(finalRoot,'_sharedWindowedData','code')); % getWindowed
dataRoot  = fullfile(finalRoot,'_sharedWindowedData','data','trained');
assert(isfolder(dataRoot), 'Trained cache not found:\n  %s', dataRoot);

outDir = fullfile(thisDir,'..','cache'); if ~exist(outDir,'dir'), mkdir(outDir); end
outMat = fullfile(outDir, sprintf('NativeAxisVsLick_%s_snap%.2fs_bin%d-%d.mat', ...
                  scheme, snapshotTime, binSizeTrials, stepSizeTrials));

%% ---------------- MAIN LOOP ----------------
files = [];
for d = 1:6
    files = [files; dir(fullfile(dataRoot, sprintf('D%d',d), '*_win.mat'))]; %#ok<AGROW>
end
assert(~isempty(files), 'No trained *_win.mat under %s', dataRoot);

TallRawTrials = table();
TsumRaw       = table();
nUsed = 0;

for fi = 1:numel(files)
    [X, info] = getWindowed(fullfile(files(fi).folder, files(fi).name), scheme);

    tc = info.tvec_centre(:);
    [~, wk] = min(abs(tc - snapshotTime));
    Xs = squeeze(X(:, :, wk));            % [nTrials x nCells] native snapshot

    lab  = double(info.Visual_Stim(:));
    resp = double(info.Response(:));
    trid = double(info.trialID(:));
    experimentID = char(string(info.session));
    behavCond    = char(string(info.day));
    animal = regexp(experimentID, '^[A-Za-z]+\d+', 'match', 'once');

    stimVals = unique(lab); stimVals = stimVals(isfinite(stimVals));
    if numel(stimVals) ~= 2 || ~any(stimVals == goStim), continue; end
    goVal   = goStim;
    nogoVal = stimVals(stimVals ~= goStim); nogoVal = nogoVal(1);

    % Go->NoGo centroid axis in native space
    goC   = mean(Xs(lab==goVal,  :), 1);
    nogoC = mean(Xs(lab==nogoVal,:), 1);
    u = (nogoC - goC); nu = norm(u); if nu < eps, continue; end
    u = u / nu;
    axisProj = (Xs - goC) * u';           % per trial

    Traw = makeTrialTable(animal, behavCond, experimentID, tc(wk), ...
                          trid, lab, resp, axisProj);
    TallRawTrials = [TallRawTrials; Traw]; %#ok<AGROW>

    Tsum = summarizeBySlidingWindow(Traw, goVal, nogoVal, binSizeTrials, stepSizeTrials, minTrialsPerStim);
    if isempty(Tsum), continue; end
    Tsum.DBname     = repmat("BatuhanAngDist",   height(Tsum), 1);
    Tsum.Animal     = repmat(string(animal),      height(Tsum), 1);
    Tsum.Behav_Cond = repmat(string(behavCond),   height(Tsum), 1);
    Tsum.Experiment = repmat(string(experimentID),height(Tsum), 1);
    Tsum.SnapshotT  = repmat(tc(wk),              height(Tsum), 1);
    TsumRaw = [TsumRaw; Tsum]; %#ok<AGROW>
    nUsed = nUsed + 1;
end

fprintf('Native: sessions used %d / %d  | snapshot ~%.2f s\n', nUsed, numel(files), snapshotTime);
save(outMat, 'TallRawTrials','TsumRaw','snapshotTime','scheme', ...
    'binSizeTrials','stepSizeTrials','minTrialsPerStim','goStim');
fprintf('Saved:\n  %s\n', outMat);

%% ========================= HELPERS =========================
function T = makeTrialTable(animal, behavCond, experimentID, snapT, trialIDs, labels, response, axisProj)
    n = numel(trialIDs);
    T = table();
    T.Animal      = repmat(string(animal),       n, 1);
    T.Behav_Cond  = repmat(string(behavCond),    n, 1);
    T.Experiment  = repmat(string(experimentID), n, 1);
    T.SnapshotT   = repmat(snapT,                n, 1);
    T.Trial       = trialIDs(:);
    T.Visual_Stim = labels(:);
    T.Response    = response(:);
    T.AxisProj    = axisProj(:);
    T.Metric      = axisProj(:);      % quantity binned against P(lick)
end

function Tsum = summarizeBySlidingWindow(Ttrials, goVal, nogoVal, binSizeTrials, stepSizeTrials, minTrialsPerStim)
    Tsum = table();
    ok = isfinite(Ttrials.Response) & isfinite(Ttrials.Metric);
    T = Ttrials(ok, :);
    if isempty(T), return; end
    for stim = [goVal, nogoVal]
        Ts = T(T.Visual_Stim == stim, :);
        if height(Ts) < minTrialsPerStim, continue; end
        d = Ts.Metric; [~, ord] = sort(d, 'ascend'); n = numel(d);
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
