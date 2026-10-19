%% pupilCapacitySplits.m
% Pupil analogue of the licking-control capacity pipeline (trained DB only).
%
% For each QC-passing session, computes a ROLLING capacity timecourse for three
% pupil-based trial partitions, in two candidate pupil windows:
%
%   Partitions (mirrors rollingCapacityLickQuantiles.m tertile logic + a new
%   Go/NoGo equalization control):
%     T1 : low  pupil tertile  (bottom third), tertiles taken WITHIN Go and
%          WITHIN NoGo separately so class balance is preserved
%     T3 : high pupil tertile  (top third), same within-class definition
%     EQ : Go/NoGo pupil-EQUALIZED subset (1:1 nearest-neighbour match on the
%          pupil metric with a caliper) -> tests whether a pupil difference
%          BETWEEN conditions drives the capacity
%
%   Windows (both analyzed per user decision):
%     baseline : [-1.5 0] s  (pre-stim tonic arousal; cleaner confound measure)
%     post     : [ 1   2] s  (mirrors LickAUCz_1_2; stimulus/choice-locked)
%
% Capacity logic is UNCHANGED from the lick pipeline: capacity_from_window()
% builds the balanced tuning_function and calls check_capacity, reading
% 2./capacity_results. Only the TRIAL SELECTION differs (pupil, not lick).
%
% Session set: the QC pass-list from pupilQualityScan.m (pupilQC_<DB>.csv,
% pass_baseline==1) MINUS the strict-exclusion sessions (marginal/glitchy).
%
% Output: one .mat per (session, window, partition) with capacity_rolling+tvec.
%
% Smoke-test overrides (set before running to preserve via clearvars):
%   RUN_SESSION_MAX, RUN_CENTERS_SUBSET, RUN_WINDOWS

clearvars -except RUN_SESSION_MAX RUN_CENTERS_SUBSET RUN_WINDOWS RUN_RESUME; clc; close all;

% Ensure the manifold-capacity library (check_capacity + deps) is on the path
if exist('check_capacity','file') == 0
    capLibRoot = fullfile(fileparts(fileparts(mfilename('fullpath'))),'toolbox_dnn-object-manifolds');
    if exist(capLibRoot,'dir'), addpath(genpath(capLibRoot)); end
    assert(exist('check_capacity','file')>0, 'check_capacity not found on path (%s)', capLibRoot);
end
% Ensure the CPLEX solver connector is on the path (cplexoptimset is a .p file,
% so exist() returns 6, not 2 -- accept any nonzero result).
if exist('cplexoptimset','file') == 0
    cplexDir = 'C:\Program Files\IBM\ILOG\CPLEX_Studio1210\cplex\matlab\x64_win64';
    if exist(cplexDir,'dir'), addpath(cplexDir); end
    assert(exist('cplexoptimset','file')>0, 'cplexoptimset (CPLEX) not found on path (%s)', cplexDir);
end

%% ---------------- SETTINGS ----------------
DB = 'BatuhanAngDist';
stringPrefOri = 'Pref_Orientation_spikes_2';   % trained DB

goVal   = 45;
noGoVal = containers.Map({'D1','D2','D3','D4','D5','D6'}, {135,90,75,70,65,60});

% Cell inclusion criteria (identical to capacity/lick pipeline)
cellCrit = struct('isCellThresh',0.8,'deconvCorrThresh',0.8, ...
    'spikeProbTCTrialsThresh',0.1,'roundnessThresh',0.2,'nPixMagFactorThresh',15);

% Time axis / windows (identical to lick pipeline)
dt        = 0.0646;
rasterLen = 140;
xval      = linspace(0,(rasterLen-1)*dt,rasterLen) - 1.5;   % 0 = stim onset
tcWindowSec = [-1 4];
winFrames   = 3;
stepFrames  = 3;
baselineWindowSec = [-1.5 0];
subtractBaseline  = true;

% Pupil metric windows (name; [start end] s rel. to stim onset)
pupilWindows = { 'baseline', [-1.5 0]; ...
                 'post',     [ 1   2] };

% Pupil preprocessing + per-trial QC (must match pupilQualityScan.m)
stimOnsetSec = 1.5;
doLowpass = true; lowpassHz = 2;
qc = struct('maxNaNFracWindow',0.10,'jumpMADk',8,'maxJumpFracTrial',0.02,'minStdRaw',1e-6);

% Partition params
minTrialsPerCondPart = 8;   % min trials/class within a partition to compute capacity
eqCaliperSD          = 0.2; % Go/NoGo match caliper, in SD of pooled pupil metric
matchSeed            = 7;   % reproducible EQ matching

% Capacity file suffix produced by check_capacity (balanced, session 1)
capacitySuffix = '_capacity_balanced_s1.mat';

% Working-set source + strict exclusions
qcCsv = fullfile(fileparts(mfilename('fullpath')), 'pupilQC', sprintf('pupilQC_%s.csv', DB));
strictExclude = { 'F10_42521_'  % D5 marginal: ~80% trials lost to glitches
                  'F07_42221_'  % D3 glitchy (already QC-fail, listed for clarity)
                  'F05_42221_'  % D3 glitchy
                  'F19_42321_' }; % D4 glitchy

% Output
saveRoot = 'D:/embedding/TraceSpikes/';   % same root the lick pipeline uses
tmpRoot  = fullfile(saveRoot, 'pupilCapTmp');   % scratch for check_capacity files
if ~exist(tmpRoot,'dir'), mkdir(tmpRoot); end

% ---- Run scope (override for smoke tests) ----
if ~exist('RUN_SESSION_MAX','var'), RUN_SESSION_MAX = Inf; end
if ~exist('RUN_CENTERS_SUBSET','var'), RUN_CENTERS_SUBSET = []; end  % [] = all timepoints
if ~exist('RUN_WINDOWS','var'), RUN_WINDOWS = {'baseline','post'}; end
if ~exist('RUN_RESUME','var'), RUN_RESUME = true; end  % skip already-saved partitions

%% ---------------- DERIVED ----------------
baselineWindow = find(xval >= baselineWindowSec(1) & xval <= baselineWindowSec(2));
tcIdx    = find(xval >= tcWindowSec(1) & xval <= tcWindowSec(2));
startIdx = tcIdx(1); endIdx = tcIdx(end);
centers  = startIdx:stepFrames:endIdx;
if ~isempty(RUN_CENTERS_SUBSET), centers = centers(RUN_CENTERS_SUBSET); end
tvec     = xval(centers);
rollingWindows = arrayfun(@(t) max(startIdx, t-winFrames+1):t, centers, 'uni', 0);

%% ---------------- WORKING SET ----------------
QCt = readtable(qcCsv);
pass = QCt(QCt.pass_baseline==1, :);
keep = ~ismember(string(pass.Experiment), string(strictExclude));
work = pass(keep, :);
fprintf('Working set: %d sessions (from %d QC-pass; %d strict-excluded).\n', ...
    height(work), height(pass), sum(~keep));

%% ---------------- CONNECT ----------------
conn = DBconnection({DB});

summary = {};
nDone = 0;
for r = 1:height(work)
    if nDone >= RUN_SESSION_MAX, break; end
    sess = char(work.Experiment(r));
    day  = char(work.BehavCond(r));
    nogoThisDay = noGoVal(day);

    fprintf('\n=== %s | %s (NoGo=%g) ===\n', day, sess, nogoThisDay);

    % Fast session-level resume: if every expected partition file already exists,
    % skip WITHOUT loading neural/pupil data (the DB query + trace decode is the
    % slow part). Only sessions needing work pay the data-loading cost.
    if RUN_RESUME
        plist = {'T1','T3','EQ'};
        sdChk = fullfile(saveRoot, ['D' day(end)], 'pupilCapacity');
        allExist = true;
        for wi = 1:numel(RUN_WINDOWS)
            for pj = 1:numel(plist)
                ff = fullfile(sdChk, sprintf('%s_%s_pupil_%s_%s_capacity_timecourse.mat', ...
                    sess, DB, RUN_WINDOWS{wi}, plist{pj}));
                if ~isfile(ff), allExist = false; break; end
            end
            if ~allExist, break; end
        end
        if allExist
            fprintf('  session fully done (%d/%d partitions present) -> skip data load\n', ...
                numel(RUN_WINDOWS)*numel(plist), numel(RUN_WINDOWS)*numel(plist));
            nDone = nDone + 1;
            continue;
        end
    end

    params = struct('DB',DB,'stringPrefOri',stringPrefOri,'cellCrit',cellCrit, ...
        'goVal',goVal,'nogoThisDay',nogoThisDay,'day',day,'sess',sess, ...
        'rasterLen',rasterLen,'endIdx',endIdx,'baselineWindow',baselineWindow, ...
        'subtractBaseline',subtractBaseline,'dt',dt,'stimOnsetSec',stimOnsetSec, ...
        'doLowpass',doLowpass,'lowpassHz',lowpassHz,'qc',qc, ...
        'pupilWindows',{pupilWindows},'runWindows',{RUN_WINDOWS}, ...
        'centers',centers,'rollingWindows',{rollingWindows},'tvec',tvec, ...
        'minTrialsPerCondPart',minTrialsPerCondPart,'eqCaliperSD',eqCaliperSD, ...
        'matchSeed',matchSeed,'capacitySuffix',capacitySuffix, ...
        'saveRoot',saveRoot,'tmpRoot',tmpRoot,'resume',RUN_RESUME);

    try
        s = runOneSession(conn, params);
        summary{end+1} = s; %#ok<AGROW>
    catch ME
        warning('Session %s failed: %s', sess, ME.message);
    end
    nDone = nDone + 1;
end

close(conn);

if ~isempty(summary)
    ST = struct2table([summary{:}]);
    outCsv = fullfile(saveRoot, sprintf('pupilCapacity_runsummary_%s.csv', DB));
    writetable(ST, outCsv);
    fprintf('\nRun summary saved: %s\n', outCsv);
    disp(ST);
end
fprintf('Done.\n');

%% ================= LOCAL FUNCTIONS =================

function out = runOneSession(conn, p)
    out = struct('Experiment',string(p.sess),'BehavCond',string(p.day), ...
        'nCells',0,'nTrialsUsedNeural',0);

    % ---- Neural traces (identical filter set to lick pipeline) ----
    nq = [' select ctt.Cell, ctt.Trial, ctt.TraceSpikes, tt.Visual_Stim ' ...
        ' from Cell_Trial_Table AS ctt ' ...
        ' inner join Trial_Table tt on tt.Experiment = ctt.Experiment AND tt.Trial=ctt.Trial' ...
        ' inner join Cell_Table ct on ct.Experiment = ctt.Experiment AND ct.Cell = ctt.Cell ' ...
        ' WHERE ct.' p.stringPrefOri ' IS NOT NULL ' ...
        ' AND (tt.Block != ''Orientation Tuning'') ' ...
        ' AND ct.isCell > '        num2str(p.cellCrit.isCellThresh) ...
        ' AND ct.DeconvCorr > '    num2str(p.cellCrit.deconvCorrThresh) ...
        ' AND ct.spikeProb_TCTrials > ' num2str(p.cellCrit.spikeProbTCTrialsThresh) ...
        ' AND ct.roundness > '     num2str(p.cellCrit.roundnessThresh) ...
        ' AND ct.nPix*ct.magFactor > ' num2str(p.cellCrit.nPixMagFactorThresh) ...
        ' AND tt.Behav_Cond = ''' p.day '''' ...
        ' AND ct.Experiment = ''' p.sess ''''];
    taskData = sortrows(select(conn, nq), {'Trial','Cell'});
    if isempty(taskData), warning('Empty neural data %s', p.sess); return; end

    uTrials = unique(taskData.Trial);
    uCells  = unique(taskData.Cell);
    nTr = numel(uTrials); nC = numel(uCells);

    T = NaN(nTr, nC, p.rasterLen, 'double');
    Vstim = NaN(nTr,1);
    for i = 1:nTr
        rows = taskData.Trial == uTrials(i);
        tc   = taskData.Cell(rows);
        raw  = taskData.TraceSpikes(rows);
        Vstim(i) = unique(taskData.Visual_Stim(rows));
        [~,loc] = ismember(tc, uCells);
        for j = 1:numel(loc)
            tr = getArrayFromByteStream(typecast(raw{j},'uint8'));
            tr = double(tr(:)'); tr = tr(1:min(end,p.rasterLen));
            if numel(tr) < p.endIdx, continue; end
            if p.subtractBaseline, tr = tr - mean(tr(p.baselineWindow),'omitnan'); end
            T(i,loc(j),1:numel(tr)) = tr;
        end
    end

    % ---- Pupil traces for the SAME session, keyed by Trial ----
    pq = sprintf([ ...
        "SELECT Trial, Pupil_Area_trace FROM Trial_Table " + ...
        "WHERE Experiment = '%s' AND (Block != 'Orientation Tuning');"], p.sess);
    Tp = fetch(conn, pq);
    pupTraces = decodePupilColumn(Tp.Pupil_Area_trace);
    pupTrialId = double(Tp.Trial);

    % Map pupil metric onto the neural trial order (uTrials)
    isGo = (Vstim == p.goVal);
    out.nCells = nC; out.nTrialsUsedNeural = nTr;

    % ---- Loop pupil windows ----
    for w = 1:size(p.pupilWindows,1)
        wname = p.pupilWindows{w,1};
        if ~any(strcmp(wname, p.runWindows)), continue; end
        wsec  = p.pupilWindows{w,2};

        % per-trial pupil metric + usable flag, aligned to uTrials
        [metric, usable] = pupilMetricForTrials(pupTraces, pupTrialId, uTrials, ...
            p.dt, p.stimOnsetSec, wsec, p.doLowpass, p.lowpassHz, p.qc);

        % z-score metric across usable trials (rank-invariant; makes caliper SD-based)
        mu = mean(metric(usable),'omitnan'); sd = std(metric(usable),'omitnan');
        if sd > 0, z = (metric - mu)/sd; else, z = metric - mu; end

        % ----- Partition trial-selection (row indices into uTrials) -----
        [T1,T3] = withinClassTertiles(z, usable, isGo);
        EQ      = matchGoNoGo(z, usable, isGo, p.eqCaliperSD, p.matchSeed);

        partSets = struct('T1',T1,'T3',T3,'EQ',EQ);
        pnames = {'T1','T3','EQ'};

        saveDir = fullfile(p.saveRoot, ['D' p.day(end)], 'pupilCapacity');
        if ~exist(saveDir,'dir'), mkdir(saveDir); end

        for pi = 1:numel(pnames)
            pn = pnames{pi};
            sel = partSets.(pn);
            nGo = sum(isGo(sel)); nNo = sum(~isGo(sel));
            out.(sprintf('n_%s_%s_Go', wname, pn)) = nGo;
            out.(sprintf('n_%s_%s_NoGo', wname, pn)) = nNo;

            if nGo < p.minTrialsPerCondPart || nNo < p.minTrialsPerCondPart
                warning('%s | %s/%s: too few trials (Go=%d NoGo=%d). Skipping.', ...
                    p.sess, wname, pn, nGo, nNo);
                continue;
            end

            outFile = fullfile(saveDir, sprintf('%s_%s_pupil_%s_%s_capacity_timecourse.mat', ...
                p.sess, p.DB, wname, pn));
            % Resumability: skip partitions already computed (unless forced)
            if p.resume && isfile(outFile)
                fprintf('  resume: skip existing %s/%s\n', wname, pn);
                continue;
            end

            Tsel = T(sel,:,:); Vsel = Vstim(sel);
            cap = NaN(numel(p.centers),1);
            for k = 1:numel(p.centers)
                prefix = fullfile(p.tmpRoot, sprintf('%s_%s_%s_%s_t%03d', ...
                    p.sess, p.DB, wname, pn, p.centers(k)));
                cap(k) = capacity_from_window(Tsel, Vsel, p.rollingWindows{k}, ...
                    prefix, p.capacitySuffix, p.goVal);
            end

            capacity_rolling = cap; tvec = p.tvec; %#ok<NASGU>
            meta = struct('DB',p.DB,'session',p.sess,'BehavCond',p.day, ...
                'pupilWindow',wname,'pupilWindowSec',wsec,'partition',pn, ...
                'nGo',nGo,'nNoGo',nNo,'eqCaliperSD',p.eqCaliperSD);
            save(outFile,'capacity_rolling','tvec','meta','-v7.3');
            fprintf('  saved %s/%s (Go=%d NoGo=%d)\n', wname, pn, nGo, nNo);
        end
    end
end

function [metric, usable] = pupilMetricForTrials(traces, trialId, uTrials, ...
    dt, stimOnsetSec, winSec, doLowpass, lowpassHz, qc)
% Returns pupil window-mean + usability, aligned to uTrials order.
    nU = numel(uTrials);
    metric = NaN(nU,1); usable = false(nU,1);
    fs = 1/dt;
    for a = 1:nU
        j = find(trialId == uTrials(a), 1);
        if isempty(j), continue; end
        v = traces{j};
        if isempty(v) || all(isnan(v)), continue; end
        v = double(v(:)'); L = numel(v);
        tRel = (0:L-1)*dt - stimOnsetSec;
        idx = tRel >= winSec(1) & tRel <= winSec(2);
        if ~any(idx), continue; end

        vm = ~isnan(v);
        if mean(vm(idx)) < (1 - qc.maxNaNFracWindow), continue; end   % coverage

        vi = v;
        if any(isnan(vi))
            g = find(~isnan(vi));
            if numel(g) < 3, continue; end
            vi = interp1(g, vi(g), 1:L, 'linear','extrap');
        end
        if std(vi(idx),0,'omitnan') <= qc.minStdRaw, continue; end     % flat

        d = diff(vi);
        if numel(d) >= 5
            m = mad(d,1);
            jf = (m>0) * mean(abs(d-median(d)) > qc.jumpMADk*m);
            if jf > qc.maxJumpFracTrial, continue; end                 % glitchy
        end

        vLP = vi;
        if doLowpass && isfinite(lowpassHz) && lowpassHz>0
            Wn = lowpassHz/(fs/2);
            if Wn < 1, [b,a2]=butter(2,Wn,'low'); vLP = filtfilt(b,a2,vi); end
        end
        metric(a) = mean(vLP(idx),'omitnan');
        usable(a) = true;
    end
end

function [T1,T3] = withinClassTertiles(z, usable, isGo)
% Bottom/top tertile taken separately within Go and within NoGo, then unioned.
    T1 = []; T3 = [];
    for cls = [true false]
        m = usable & (isGo==cls);
        idx = find(m);
        if numel(idx) < 3, continue; end
        vals = z(idx);
        q33 = quantile(vals, 1/3); q67 = quantile(vals, 2/3);
        T1 = [T1; idx(vals <= q33)]; %#ok<AGROW>
        T3 = [T3; idx(vals >= q67)]; %#ok<AGROW>
    end
    T1 = sort(T1); T3 = sort(T3);
end

function EQ = matchGoNoGo(z, usable, isGo, caliperSD, seed)
% 1:1 nearest-neighbour match (without replacement) between Go and NoGo on the
% pupil metric z, keeping pairs within caliperSD. Returns unioned row indices.
    rng(seed);
    goIdx = find(usable &  isGo);
    noIdx = find(usable & ~isGo);
    EQ = [];
    if isempty(goIdx) || isempty(noIdx), return; end
    zg = z(goIdx); zn = z(noIdx);
    usedNo = false(numel(noIdx),1);
    order = randperm(numel(goIdx));
    for ii = order
        d = abs(zn - zg(ii));
        d(usedNo) = Inf;
        [dm, k] = min(d);
        if isfinite(dm) && dm <= caliperSD
            usedNo(k) = true;
            EQ = [EQ; goIdx(ii); noIdx(k)]; %#ok<AGROW>
        end
    end
    EQ = sort(EQ);
end

function traces = decodePupilColumn(rawCol)
    n = numel(rawCol); traces = cell(n,1);
    if ~iscell(rawCol), rawCol = num2cell(rawCol); end
    for i = 1:n
        raw = rawCol{i};
        if iscell(raw) && numel(raw)==1, raw = raw{1}; end
        if isempty(raw), traces{i} = []; continue; end
        try
            tr = getArrayFromByteStream(typecast(raw,'uint8'));
            traces{i} = double(tr(:)');
        catch
            traces{i} = [];
        end
    end
end

function cap = capacity_from_window(T, Vstim, winIdx, name_prefix, capacity_suffix, goVal)
% VERBATIM capacity logic from rollingCapacityLickQuantiles.m (only the Go value
% is passed in, and stale output files are cleaned to avoid check_capacity's
% skip-existing behavior returning a cached value).
    cap = NaN;
    if isempty(winIdx), return; end
    X = squeeze(nansum(T(:,:,winIdx),3));        % [nTrials x nCells]
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
    data_titles     = {name_prefix}; %#ok<NASGU>
    session_ids     = 1;

    % Clean any stale files so check_capacity recomputes (it skips if present)
    tuneFile = [name_prefix '_tuning.mat'];
    capFile  = [name_prefix capacity_suffix];
    if isfile(tuneFile), delete(tuneFile); end
    if isfile(capFile),  delete(capFile);  end

    save(tuneFile, 'tuning_function', 'data_titles', '-v7.3');
    try
        check_capacity(name_prefix, session_ids);
    catch ME
        warning('check_capacity failed at %s: %s', name_prefix, ME.message);
        return;
    end
    if ~isfile(capFile)
        d = dir([fileparts(name_prefix) filesep '*' capacity_suffix]);
        if isempty(d), return; else, capFile = fullfile(d(1).folder, d(1).name); end
    end
    file = matfile(capFile);
    cap = 2 ./ file.capacity_results;
end
