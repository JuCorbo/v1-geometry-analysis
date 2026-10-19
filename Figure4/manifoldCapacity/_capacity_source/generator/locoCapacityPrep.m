%% locoCapacityPrep.m  — PHASE 1 of 2 (data prep; DB, no capacity)
% Builds and SAVES, per working-set session, everything the capacity step needs:
%   T (decoded, baseline-subtracted neural [nTr x nC x rasterLen]), Vstim,
%   the T1/T3/EQ trial-index sets (post [1 2]s locomotion metric), and the
%   rolling-window schedule. NO check_capacity here.
%
% WHY TWO PHASES: check_capacity (CPLEX MEX) segfaults when called inside the
% monolithic fetch+compute pipeline, but runs reliably in a clean process that
% only loads prepared arrays, so Phase 2 (locoCapacityCompute.m) runs
% check_capacity in a fresh, minimal process.
%
% Output: D:/embedding/TraceSpikes/locoPrep/<sess>_<DB>_prep.mat
%
% Smoke override: RUN_SESSION_MAX

clearvars -except RUN_SESSION_MAX; clc; close all;

%% ---------------- SETTINGS (must match locoCapacityCompute.m) ----------------
DB = 'BatuhanAngDist';
stringPrefOri = 'Pref_Orientation_spikes_2';
goVal   = 45;
noGoVal = containers.Map({'D1','D2','D3','D4','D5','D6'}, {135,90,75,70,65,60});
cellCrit = struct('isCellThresh',0.8,'deconvCorrThresh',0.8, ...
    'spikeProbTCTrialsThresh',0.1,'roundnessThresh',0.2,'nPixMagFactorThresh',15);

dt        = 0.0646;
rasterLen = 140;
xval      = linspace(0,(rasterLen-1)*dt,rasterLen) - 1.5;
tcWindowSec = [-1 4];
winFrames   = 3;
stepFrames  = 3;
baselineWindowSec = [-1.5 0];
subtractBaseline  = true;

locoWinSec   = [1 2];          % single post window
stimOnsetSec = 1.5;
locoLowpassHz = 1;
locoMinFrames = 60;
eqCaliperSD  = 0.2;
matchSeed    = 7;

qcCsv           = fullfile(fileparts(mfilename('fullpath')), 'locoQC', sprintf('locoQC_%s.csv', DB));
minSessStdSplit = 3;

prepRoot = 'D:/embedding/TraceSpikes/locoPrep';
if ~exist(prepRoot,'dir'), mkdir(prepRoot); end

if ~exist('RUN_SESSION_MAX','var'), RUN_SESSION_MAX = Inf; end

%% ---------------- DERIVED (shared schedule) ----------------
baselineWindow = find(xval >= baselineWindowSec(1) & xval <= baselineWindowSec(2));
tcIdx    = find(xval >= tcWindowSec(1) & xval <= tcWindowSec(2));
startIdx = tcIdx(1); endIdx = tcIdx(end);
centers  = startIdx:stepFrames:endIdx;
tvec     = xval(centers);
rollingWindows = arrayfun(@(t) max(startIdx, t-winFrames+1):t, centers, 'uni', 0);

%% ---------------- WORKING SET ----------------
QCt  = readtable(qcCsv);
pass = QCt(QCt.pass_post==1, :);
work = pass(pass.sessStd >= minSessStdSplit, :);
fprintf('Working set: %d sessions.\n', height(work));

conn = DBconnection({DB});
nDone = 0;
for r = 1:height(work)
    if nDone >= RUN_SESSION_MAX, break; end
    sess = char(work.Experiment(r)); day = char(work.BehavCond(r));
    nogoThisDay = noGoVal(day);
    outFile = fullfile(prepRoot, sprintf('%s_%s_prep.mat', sess, DB));
    fprintf('\n=== PREP %s | %s ===\n', day, sess);

    % ---- neural ----
    nq = [' select ctt.Cell, ctt.Trial, ctt.TraceSpikes, tt.Visual_Stim ' ...
        ' from Cell_Trial_Table AS ctt ' ...
        ' inner join Trial_Table tt on tt.Experiment = ctt.Experiment AND tt.Trial=ctt.Trial' ...
        ' inner join Cell_Table ct on ct.Experiment = ctt.Experiment AND ct.Cell = ctt.Cell ' ...
        ' WHERE ct.' stringPrefOri ' IS NOT NULL AND (tt.Block != ''Orientation Tuning'') ' ...
        ' AND ct.isCell > '        num2str(cellCrit.isCellThresh) ...
        ' AND ct.DeconvCorr > '    num2str(cellCrit.deconvCorrThresh) ...
        ' AND ct.spikeProb_TCTrials > ' num2str(cellCrit.spikeProbTCTrialsThresh) ...
        ' AND ct.roundness > '     num2str(cellCrit.roundnessThresh) ...
        ' AND ct.nPix*ct.magFactor > ' num2str(cellCrit.nPixMagFactorThresh) ...
        ' AND tt.Behav_Cond = ''' day ''' AND ct.Experiment = ''' sess ''''];
    taskData = sortrows(select(conn, nq), {'Trial','Cell'});
    if isempty(taskData), warning('Empty neural %s', sess); continue; end
    uTrials = unique(taskData.Trial); uCells = unique(taskData.Cell);
    nTr = numel(uTrials); nC = numel(uCells);
    T = NaN(nTr, nC, rasterLen);
    Vstim = NaN(nTr,1);
    for i = 1:nTr
        rows = taskData.Trial == uTrials(i);
        tc = taskData.Cell(rows); raw = taskData.TraceSpikes(rows);
        Vstim(i) = unique(taskData.Visual_Stim(rows));
        [~,loc] = ismember(tc, uCells);
        for j = 1:numel(loc)
            tr = getArrayFromByteStream(typecast(raw{j},'uint8'));
            tr = double(tr(:)'); tr = tr(1:min(end,rasterLen));
            if numel(tr) < endIdx, continue; end
            if subtractBaseline, tr = tr - mean(tr(baselineWindow),'omitnan'); end
            T(i,loc(j),1:numel(tr)) = tr;
        end
    end

    % ---- locomotion metric (post [1 2]s) ----
    lq = sprintf(["SELECT Trial, Locomotion_trace FROM Trial_Table " + ...
        "WHERE Experiment = '%s' AND (Block != 'Orientation Tuning');"], sess);
    Tl = fetch(conn, lq);
    locoTraces = decodeLocoColumn(Tl.Locomotion_trace);
    locoTrialId = double(Tl.Trial);
    isGo = (Vstim == goVal);
    [metric, usable] = locoMetricForTrials(locoTraces, locoTrialId, uTrials, ...
        dt, stimOnsetSec, locoWinSec, locoLowpassHz, locoMinFrames);
    mu = mean(metric(usable),'omitnan'); sd = std(metric(usable),'omitnan');
    if sd > 0, z = (metric - mu)/sd; else, z = metric - mu; end

    [T1,T3] = withinClassTertiles(z, usable, isGo);
    EQ      = matchGoNoGo(z, usable, isGo, eqCaliperSD, matchSeed);

    meta = struct('DB',DB,'session',sess,'BehavCond',day,'goVal',goVal, ...
        'nogoThisDay',nogoThisDay,'locoWinSec',locoWinSec,'eqCaliperSD',eqCaliperSD);
    save(outFile, 'T','Vstim','isGo','z','usable','T1','T3','EQ', ...
        'centers','rollingWindows','tvec','goVal','nogoThisDay','day','sess','DB','meta','-v7.3');
    fprintf('  saved prep: nTr=%d nC=%d | T1=%d T3=%d EQ=%d -> %s\n', ...
        nTr, nC, numel(T1), numel(T3), numel(EQ), outFile);
    nDone = nDone + 1;
end
close(conn);
fprintf('\nPrep done (%d sessions).\n', nDone);

%% ================= LOCAL FUNCTIONS =================
function traces = decodeLocoColumn(rawCol)
    n = numel(rawCol); traces = cell(n,1);
    if ~iscell(rawCol), rawCol = num2cell(rawCol); end
    for i = 1:n
        raw = rawCol{i};
        if iscell(raw) && numel(raw)==1, raw = raw{1}; end
        if isempty(raw), traces{i} = []; continue; end
        try, tr = getArrayFromByteStream(typecast(raw,'uint8')); traces{i} = double(tr(:)');
        catch, traces{i} = []; end
    end
end
function [metric, usable] = locoMetricForTrials(traces, trialId, uTrials, dt, stimOnsetSec, winSec, lowpassHz, minFrames)
    nU = numel(uTrials); metric = NaN(nU,1); usable = false(nU,1); fs = 1/dt;
    for a = 1:nU
        j = find(trialId == uTrials(a), 1); if isempty(j), continue; end
        v = traces{j}; if isempty(v), continue; end
        v = double(v(:)'); L = numel(v); if L < minFrames, continue; end
        v(~isfinite(v)) = 0;
        tRel = (0:L-1)*dt - stimOnsetSec; idx = tRel >= winSec(1) & tRel <= winSec(2);
        if ~any(idx), continue; end
        if L > 12, vLP = lowpass(v, lowpassHz, fs); else, vLP = v; end
        metric(a) = mean(vLP(idx),'omitnan'); usable(a) = true;
    end
end
function [T1,T3] = withinClassTertiles(z, usable, isGo)
    T1 = []; T3 = [];
    for cls = [true false]
        idx = find(usable & (isGo==cls)); if numel(idx) < 3, continue; end
        vals = z(idx); q33 = quantile(vals,1/3); q67 = quantile(vals,2/3);
        T1 = [T1; idx(vals <= q33)]; T3 = [T3; idx(vals >= q67)]; %#ok<AGROW>
    end
    T1 = sort(T1); T3 = sort(T3);
end
function EQ = matchGoNoGo(z, usable, isGo, caliperSD, seed)
    rng(seed); goIdx = find(usable & isGo); noIdx = find(usable & ~isGo); EQ = [];
    if isempty(goIdx) || isempty(noIdx), return; end
    zg = z(goIdx); zn = z(noIdx); usedNo = false(numel(noIdx),1); order = randperm(numel(goIdx));
    for ii = order
        d = abs(zn - zg(ii)); d(usedNo) = Inf; [dm, k] = min(d);
        if isfinite(dm) && dm <= caliperSD, usedNo(k) = true; EQ = [EQ; goIdx(ii); noIdx(k)]; end %#ok<AGROW>
    end
    EQ = sort(EQ);
end
