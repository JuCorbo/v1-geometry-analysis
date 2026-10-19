%% Super-mouse capacity timecourse (TRAINED only) with trial sorting by LICK vs NO-LICK
% UPDATED per your requests (FINALIZED):
%   - Trial sorting window declared at top (default = [2 3] s).
%     If sortWinSec == [2 3] and useResponseIfMatches=true, use Trial_Table.Response.
%     Otherwise compute lick-in-window from Lick_trace/Lick_dt with per-session bimodal threshold.
%   - Neuron query includes: ct.Pref_Orientation_spikes_2 IS NOT NULL
%   - NO averaged pseudotrials.
%   - CRITICAL UPDATE: pseudo-trial trial IDs are sampled ONCE per replicate+mode
%     and reused across ALL timepoints (so you “use all timepoints for the selected trials”).
%   - Keeps SAME NEURONS across modes (all / lick / nolick) within a replicate.
%   - Saves pseudo-trial lick traces (weighted by #neurons contributing from each session).
%   - Uses ONLY rolling windows (expanding removed).

clear; clc; close all;

%% -------------------- USER SETTINGS --------------------
rootDir    = 'D:/embedding/TraceSpikes/';
dbTrained  = 'BatuhanAngDist';

BehavCond  = {'D1','D2','D3','D4','D5','D6'};
goVal      = 45;
noGoVal    = [135, 90, 75, 70, 65, 60];

% ---------- Trial sorting window (relative to stim onset) ----------
sortWinSec = [1 2];     % DEFAULT = [2 3] seconds (matches Trial_Table.Response definition)

% If sortWinSec == [2 3], use Response as lick/no-lick flag (fast and exact).
useResponseIfMatches = true;

% If sortWinSec differs from [2 3], compute lick presence from Lick_trace/Lick_dt
stimOnsetSec  = 1.5;     % alignment for lick trace: tRel = tAbs - stimOnsetSec
doSmoothLick   = true;   % smoothing before thresholding
smoothWinSec   = 0.15;   % smoothing window (seconds)
lickMinDurSec  = 0.00;   % require >= this contiguous duration above threshold (0 = any sample)

% Per-session threshold from bimodal lick distribution
threshClipPrct = [1 99];     % robust clip before fitting
threshMaxPts   = 30000;      % downsample points for threshold fit

% ---------- Supermouse resampling ----------
nPseudoSessions       = 5;     % per day
nNeuronsTarget        = 500;   % neurons per pseudo-session
nPseudoTrialsPerCond  = 100;   % number of RAW pseudo-trials per cond (Go/NoGo)
minTrialsPerCond      = 8;     % min eligible REAL trials per (session,cond,mode) to allow drawing

% If you want to require that EVERY contributing session has >=minTrialsPerCond,
% set strictPerSessionEligibility=true. If false, sessions with low eligible counts
% are still allowed; missing draws become NaN and those neurons are left as NaN.
strictPerSessionEligibility = false;

% ---------- Time axis / windows ----------
dt        = 0.0646;                 % s per frame
rasterLen = 140;
xval      = linspace(0,(rasterLen-1)*dt,rasterLen) - 1.5;  % 0 = stim onset

tcWindowSec = [-1 +4];
winFrames   = 3;
stepFrames  = 3;

% Baseline subtraction (spikes)
baselineWindowSec = [-1.5 0];
subtractBaseline  = true;

% Cell inclusion criteria (your standard)
cellCrit = struct();
cellCrit.isCellThresh            = 0.8;
cellCrit.deconvCorrThresh        = 0.8;
cellCrit.spikeProbTCTrialsThresh = 0.1;
cellCrit.roundnessThresh         = 0.2;
cellCrit.nPixMagFactorThresh     = 15;

% Require tuning existence (requested)
requirePrefOri2 = true; % adds: ct.Pref_Orientation_spikes_2 IS NOT NULL

% Capacity file suffix produced by your toolbox
capacitySuffix = '_capacity_balanced_s1.mat';

% Plot colors (for per-day summary)
colAllRoll     = [0.20 0.20 0.20];
colLickRoll    = [0.90 0.20 0.20];
colNoLickRoll  = [0.20 0.40 0.90];

% Save outputs
saveOutputs = true;

%% -------------------- DERIVED INDICES --------------------
baselineWindow = find(xval >= baselineWindowSec(1) & xval <= baselineWindowSec(2));
tcIdx     = find(xval >= tcWindowSec(1) & xval <= tcWindowSec(2));
startIdx  = tcIdx(1);
endIdx    = tcIdx(end);
centers   = startIdx:stepFrames:endIdx;
tvec      = xval(centers);

rollingWindows = arrayfun(@(t) max(startIdx, t-winFrames+1):t, centers, 'uni', 0);

%% -------------------- DAY LOOP --------------------
dayFolders = dir(fullfile(rootDir, 'D*'));
dayFolders = dayFolders([dayFolders.isdir]);
if isempty(dayFolders), error('No day folders (D*) found under %s', rootDir); end

for d = 1:numel(dayFolders)
    day = dayFolders(d).name;           % e.g., 'D1'
    dayNum = sscanf(day, 'D%d');
    if isempty(dayNum) || dayNum < 1 || dayNum > 6, continue; end

    nogoThisDay = noGoVal(dayNum);

    fprintf('\n\n========================= %s (NoGo=%g) =========================\n', day, nogoThisDay);

    outRoot = fullfile(rootDir, day, 'supermouse');
    if saveOutputs && ~exist(outRoot,'dir'), mkdir(outRoot); end

    % Pool session data for this day (TRAINED only)
    [grp, nCellsTotal, sortMeta] = load_group_pool_trained_sorted( ...
        rootDir, day, dbTrained, rasterLen, subtractBaseline, baselineWindow, endIdx, cellCrit, requirePrefOri2, ...
        sortWinSec, useResponseIfMatches, stimOnsetSec, doSmoothLick, smoothWinSec, lickMinDurSec, ...
        threshClipPrct, threshMaxPts);

    if isempty(grp.sessionsT)
        warning('%s: no sessions pooled. Skipping.', day);
        continue;
    end

    nSess = numel(grp.sessionsT);
    fprintf('%s: pooled %d sessions, totalCells=%d\n', day, nSess, nCellsTotal);
    fprintf('Sorting window = [%g %g] s, method = %s\n', sortWinSec(1), sortWinSec(2), sortMeta.method);

    % Precompute lick grids per session once (for weighted pseudo-lick traces)
    [sessionsLickGrid, tLickGrid] = build_sessionsLickGrid(grp, stimOnsetSec, doSmoothLick, smoothWinSec);

    % Prepare outputs: rolling capacity only
    capAll    = NaN(nPseudoSessions, numel(centers));
    capLick   = NaN(nPseudoSessions, numel(centers));
    capNoLick = NaN(nPseudoSessions, numel(centers));

    auxAll    = cell(nPseudoSessions,1);
    auxLick   = cell(nPseudoSessions,1);
    auxNoLick = cell(nPseudoSessions,1);

    % Session cell counts for neuron sampling
    sessCellCounts = cellfun(@(Ts) size(Ts,2), grp.sessionsT);

    for r = 1:nPseudoSessions
        fprintf('  Pseudo-session %d/%d...\n', r, nPseudoSessions);

        % 1) sample neuron identities across sessions ONCE
        sel = sample_neurons_across_sessions(sessCellCounts, nNeuronsTarget);

        % Build per-session neuron column mapping
        [sessList, colsBySess, nNeurBySess, nNeurTotal, wBySess] = build_session_column_map(sel, nSess);

        % 2) draw trial IDs ONCE per mode (reused across all timepoints)
        drawAll_r    = draw_trial_ids_once(grp.sessionsV, grp.sessionsHasLick, sessList, goVal, nogoThisDay, 'all',    nPseudoTrialsPerCond, minTrialsPerCond, strictPerSessionEligibility);
        drawLick_r   = draw_trial_ids_once(grp.sessionsV, grp.sessionsHasLick, sessList, goVal, nogoThisDay, 'lick',   nPseudoTrialsPerCond, minTrialsPerCond, strictPerSessionEligibility);
        drawNoLick_r = draw_trial_ids_once(grp.sessionsV, grp.sessionsHasLick, sessList, goVal, nogoThisDay, 'nolick', nPseudoTrialsPerCond, minTrialsPerCond, strictPerSessionEligibility);

        % Build weighted pseudo-lick traces ONCE per mode (fixed across time)
        [lickPseudoAllGo,   lickPseudoAllNoGo]    = build_weighted_pseudo_lick(sessionsLickGrid, drawAll_r,    wBySess, tLickGrid);
        [lickPseudoLickGo,  lickPseudoLickNoGo]   = build_weighted_pseudo_lick(sessionsLickGrid, drawLick_r,   wBySess, tLickGrid);
        [lickPseudoNoGoGo,  lickPseudoNoGoNoGo]   = build_weighted_pseudo_lick(sessionsLickGrid, drawNoLick_r, wBySess, tLickGrid);

        % Store aux (once)
        auxAll{r} = struct('sel',sel,'sessList',sessList,'colsBySess',{colsBySess},'nNeurBySess',nNeurBySess,'nNeurTotal',nNeurTotal, ...
                           'draw',drawAll_r,'tLickGrid',tLickGrid,'lickGo',lickPseudoAllGo,'lickNoGo',lickPseudoAllNoGo);
        auxLick{r}= struct('sel',sel,'sessList',sessList,'colsBySess',{colsBySess},'nNeurBySess',nNeurBySess,'nNeurTotal',nNeurTotal, ...
                           'draw',drawLick_r,'tLickGrid',tLickGrid,'lickGo',lickPseudoLickGo,'lickNoGo',lickPseudoLickNoGo);
        auxNoLick{r}=struct('sel',sel,'sessList',sessList,'colsBySess',{colsBySess},'nNeurBySess',nNeurBySess,'nNeurTotal',nNeurTotal, ...
                           'draw',drawNoLick_r,'tLickGrid',tLickGrid,'lickGo',lickPseudoNoGoGo,'lickNoGo',lickPseudoNoGoNoGo);

        % 3) time loop: ONLY Xi depends on time; trials are fixed by draw*
        for k = 1:numel(centers)
            winIdx = rollingWindows{k};

            % Build per-session Xi blocks for current window, for selected neurons
            Xi = build_Xi_blocks(grp.sessionsT, sel, winIdx);

            % Compute capacities for each mode using fixed draws
            tagAll    = sprintf('%s_All_roll_t%03d_rep%02d', day, centers(k), r);
            tagLick   = sprintf('%s_Lick_roll_t%03d_rep%02d', day, centers(k), r);
            tagNoLick = sprintf('%s_NoLick_roll_t%03d_rep%02d', day, centers(k), r);

            capAll(r,k)    = capacity_from_fixed_draw(Xi, colsBySess, sessList, goVal, nogoThisDay, drawAll_r,    outRoot, tagAll,    capacitySuffix);
            capLick(r,k)   = capacity_from_fixed_draw(Xi, colsBySess, sessList, goVal, nogoThisDay, drawLick_r,   outRoot, tagLick,   capacitySuffix);
            capNoLick(r,k) = capacity_from_fixed_draw(Xi, colsBySess, sessList, goVal, nogoThisDay, drawNoLick_r, outRoot, tagNoLick, capacitySuffix);
        end
    end

    % ---- Plot summary per day ----
    fig = figure('Color','w','Position',[120 120 1100 620]); hold on;
    plot_group_mean_sem(tvec, capAll,    colAllRoll,    'All');
    plot_group_mean_sem(tvec, capLick,   colLickRoll,   'Lick');
    plot_group_mean_sem(tvec, capNoLick, colNoLickRoll, 'No-lick');
    xline(0,'--k','Stim onset','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
    xlabel('Time from stimulus onset (s)');
    ylabel('Capacity (2 / capacity\_results)');
    title(sprintf(['Super-mouse TRAINED — %s — NoGo=%g | %d neurons, %d pseudo-trials/cond, %d pseudo-sessions\n' ...
                   'Trial IDs fixed across time | Trial sorting window [%g %g] s (%s)'], ...
        day, nogoThisDay, nNeuronsTarget, nPseudoTrialsPerCond, nPseudoSessions, sortWinSec(1), sortWinSec(2), sortMeta.method), ...
        'Interpreter','none');
    legend('Location','best'); grid on; box off;

    if saveOutputs
        outPng = fullfile(outRoot, sprintf('supermouse_%s_trained_sortWin%_FIXEDTRIALS_summary.png', day, num2str(sortWinSec)));
        outMat = fullfile(outRoot, sprintf('supermouse_%s_trained_sortWin%_FIXEDTRIALS_summary.mat', day,num2str(sortWinSec)));
        saveas(fig, outPng);

        save(outMat, ...
            'capAll','capLick','capNoLick','auxAll','auxLick','auxNoLick', ...
            'tvec','day','dbTrained','goVal','nogoThisDay', ...
            'nNeuronsTarget','nPseudoTrialsPerCond','nPseudoSessions', ...
            'winFrames','stepFrames','tcWindowSec','baselineWindowSec','subtractBaseline','cellCrit','requirePrefOri2', ...
            'sortWinSec','useResponseIfMatches','stimOnsetSec','doSmoothLick','smoothWinSec','lickMinDurSec', ...
            'threshClipPrct','threshMaxPts','sortMeta','strictPerSessionEligibility');
        fprintf('Saved outputs:\n  %s\n  %s\n', outPng, outMat);
    end

    close(fig);
end

disp('All days processed.');

%% ================== FUNCTIONS ==================

function [grp, totalCells, sortMeta] = load_group_pool_trained_sorted( ...
    rootDir, day, dbName, rasterLen, subtractBaseline, baselineWindow, endIdx, cellCrit, requirePrefOri2, ...
    sortWinSec, useResponseIfMatches, stimOnsetSec, doSmoothLick, smoothWinSec, lickMinDurSec, ...
    clipPrct, maxPts)

    dayDir = fullfile(rootDir, day, 'timecourse');
    if ~exist(dayDir,'dir')
        warning('Day folder not found: %s', dayDir);
        grp = struct('sessionsT',{{}},'sessionsV',{{}},'sessionsHasLick',{{}},'sessionsTrialID',{{}},'sessionsExpName',{{}},'sessionsLickRaw',{{}},'sessionsLickDt',{{}});
        totalCells = 0;
        sortMeta = struct('method','none','note','day folder missing');
        return;
    end

    summaries = dir(fullfile(dayDir, sprintf('*_%s_capacity_timecourse.mat', dbName)));
    if isempty(summaries)
        warning('No session summaries for %s %s. Pool will be empty.', day, dbName);
        grp = struct('sessionsT',{{}},'sessionsV',{{}},'sessionsHasLick',{{}},'sessionsTrialID',{{}},'sessionsExpName',{{}},'sessionsLickRaw',{{}},'sessionsLickDt',{{}});
        totalCells = 0;
        sortMeta = struct('method','none','note','no session summaries');
        return;
    end

    % Decide sorting method
    useRespShortcut = false;
    if useResponseIfMatches
        useRespShortcut = (numel(sortWinSec)==2) && abs(sortWinSec(1)-2)<1e-9 && abs(sortWinSec(2)-3)<1e-9;
    end
    if useRespShortcut
        sortMeta = struct('method','Response (fast path)','note','Response==1 means lick in 2–3 s by DB definition');
    else
        sortMeta = struct('method','Lick_trace (computed)','note','Computed lick-in-window from analog trace + session bimodal threshold');
    end

    prefClause = '';
    if requirePrefOri2
        prefClause = ' AND ct.Pref_Orientation_spikes_2 IS NOT NULL ';
    end

    conn = DBconnection({dbName});

    sessionsT = {};
    sessionsV = {};
    sessionsHasLick = {};
    sessionsTrialID = {};
    sessionsExpName = {};
    sessionsLickRaw = {};   % per trial (cell array)
    sessionsLickDt  = {};   % per trial (vector)
    totalCells = 0;

    for k = 1:numel(summaries)
        % session id used in DB = before "_<dbName>_capacity_timecourse.mat"
        sessName = extractBefore(summaries(k).name, ['_' dbName '_capacity_timecourse.mat']);

        if useRespShortcut
            selectquery = [ ...
                ' select ctt.Experiment,' ...
                ' ctt.Cell, ctt.Trial, ctt.TraceSpikes,' ...
                ' tt.Visual_Stim, tt.Response ' ...
                ' from Cell_Trial_Table AS ctt ' ...
                ' inner join Trial_Table tt ' ...
                ' on tt.Experiment = ctt.Experiment AND tt.Trial=ctt.Trial' ...
                ' inner join Cell_Table ct ' ...
                ' on ct.Experiment = ctt.Experiment AND ct.Cell = ctt.Cell ' ...
                ' WHERE ct.isCell > ' num2str(cellCrit.isCellThresh) ...
                ' AND ct.DeconvCorr > ' num2str(cellCrit.deconvCorrThresh) ...
                ' AND ct.spikeProb_TCTrials > ' num2str(cellCrit.spikeProbTCTrialsThresh) ...
                ' AND ct.roundness > ' num2str(cellCrit.roundnessThresh) ...
                ' AND ct.nPix*ct.magFactor > ' num2str(cellCrit.nPixMagFactorThresh) ...
                prefClause ...
                ' AND (tt.Block != ''Orientation Tuning'') ' ...
                ' AND ct.Experiment= ''' sessName '''' ];
        else
            selectquery = [ ...
                ' select ctt.Experiment,' ...
                ' ctt.Cell, ctt.Trial, ctt.TraceSpikes,' ...
                ' tt.Visual_Stim, tt.Response, tt.Lick_trace, tt.Lick_dt ' ...
                ' from Cell_Trial_Table AS ctt ' ...
                ' inner join Trial_Table tt ' ...
                ' on tt.Experiment = ctt.Experiment AND tt.Trial=ctt.Trial' ...
                ' inner join Cell_Table ct ' ...
                ' on ct.Experiment = ctt.Experiment AND ct.Cell = ctt.Cell ' ...
                ' WHERE ct.isCell > ' num2str(cellCrit.isCellThresh) ...
                ' AND ct.DeconvCorr > ' num2str(cellCrit.deconvCorrThresh) ...
                ' AND ct.spikeProb_TCTrials > ' num2str(cellCrit.spikeProbTCTrialsThresh) ...
                ' AND ct.roundness > ' num2str(cellCrit.roundnessThresh) ...
                ' AND ct.nPix*ct.magFactor > ' num2str(cellCrit.nPixMagFactorThresh) ...
                prefClause ...
                ' AND (tt.Block != ''Orientation Tuning'') ' ...
                ' AND ct.Experiment= ''' sessName '''' ];
        end

        S = select(conn, selectquery);
        if isempty(S)
            warning('No task data for %s (%s)', sessName, dbName);
            continue;
        end
        S = sortrows(S, {'Trial','Cell'});

        unique_trials = unique(S.Trial);
        unique_cells  = unique(S.Cell);
        nTrials = numel(unique_trials);
        nCells  = numel(unique_cells);

        T = NaN(nTrials, nCells, rasterLen, 'double');
        Vstim = NaN(nTrials,1);
        Resp  = NaN(nTrials,1);

        lickTraceByTrial = cell(nTrials,1);
        lickDtByTrial    = NaN(nTrials,1);

        for i = 1:nTrials
            rows = S.Trial == unique_trials(i);
            trial_cells = S.Cell(rows);
            rawTraces   = S.TraceSpikes(rows);

            thisStim = unique(S.Visual_Stim(rows));
            thisResp = unique(S.Response(rows));
            if numel(thisStim) ~= 1 || numel(thisResp) ~= 1
                warning('Non-unique Visual_Stim/Response within trial %d in %s. Skipping trial.', unique_trials(i), sessName);
                continue;
            end
            Vstim(i) = double(thisStim);
            Resp(i)  = double(thisResp);

            if ~useRespShortcut
                rawL  = unique(S.Lick_trace(rows));
                rawDt = unique(S.Lick_dt(rows));
                if numel(rawL)==1 && numel(rawDt)==1 && ~isempty(rawL{1})
                    lickTraceByTrial{i} = decodeOneByteStreamTrace(rawL{1});
                    lickDtByTrial(i) = double(rawDt);
                else
                    lickTraceByTrial{i} = [];
                    lickDtByTrial(i) = NaN;
                end
            end

            [~, loc] = ismember(trial_cells, unique_cells);
            for j = 1:numel(loc)
                tr = getArrayFromByteStream(typecast(rawTraces{j}, 'uint8'));
                tr = double(tr(:)');
                tr = tr(1:min(end, rasterLen));
                if numel(tr) < endIdx, continue; end
                if subtractBaseline
                    tr = tr - mean(tr(baselineWindow), 'omitnan');
                end
                T(i, loc(j), 1:numel(tr)) = tr;
            end
        end

        okT = isfinite(Vstim) & isfinite(Resp);
        T = T(okT,:,:);
        Vstim = Vstim(okT);
        Resp  = Resp(okT);
        unique_trials = unique_trials(okT);

        if ~useRespShortcut
            lickTraceByTrial = lickTraceByTrial(okT);
            lickDtByTrial    = lickDtByTrial(okT);
        else
            lickTraceByTrial = cell(sum(okT),1);
            lickDtByTrial    = NaN(sum(okT),1);
        end

        if isempty(T)
            warning('Session %s had no valid trials after filtering.', sessName);
            continue;
        end

        if useRespShortcut
            hasLick = (Resp == 1);
        else
            lickDtSess = median(lickDtByTrial(isfinite(lickDtByTrial)));
            if ~(isfinite(lickDtSess) && lickDtSess > 0)
                warning('Session %s: invalid Lick_dt; cannot compute lick window. Skipping session.', sessName);
                continue;
            end

            lickThresh = computeBimodalThreshold_session(lickTraceByTrial, doSmoothLick, smoothWinSec, lickDtSess, clipPrct, maxPts);
            if ~isfinite(lickThresh)
                warning('Session %s: failed to compute lick threshold. Skipping session.', sessName);
                continue;
            end

            hasLick = computeHasLickInWindow(lickTraceByTrial, lickDtSess, stimOnsetSec, sortWinSec, doSmoothLick, smoothWinSec, lickThresh, lickMinDurSec);
        end

        sessionsT{end+1} = T; %#ok<AGROW>
        sessionsV{end+1} = Vstim; %#ok<AGROW>
        sessionsHasLick{end+1} = hasLick; %#ok<AGROW>
        sessionsTrialID{end+1} = unique_trials; %#ok<AGROW>
        sessionsExpName{end+1} = sessName; %#ok<AGROW>
        sessionsLickRaw{end+1} = lickTraceByTrial; %#ok<AGROW>
        sessionsLickDt{end+1}  = lickDtByTrial; %#ok<AGROW>

        totalCells = totalCells + nCells;
    end

    close(conn);

    grp = struct('sessionsT', {sessionsT}, 'sessionsV', {sessionsV}, 'sessionsHasLick', {sessionsHasLick}, ...
                 'sessionsTrialID',{sessionsTrialID}, 'sessionsExpName',{sessionsExpName}, ...
                 'sessionsLickRaw',{sessionsLickRaw}, 'sessionsLickDt',{sessionsLickDt});
end

function tr = decodeOneByteStreamTrace(raw)
    try
        tr = getArrayFromByteStream(typecast(raw, 'uint8'));
        tr = double(tr(:)'); % row
    catch
        tr = [];
    end
end

function thr = computeBimodalThreshold_session(traces, doSmooth, smoothWinSec, dt, clipPrct, maxPts)
% Pool samples across trials; optional smoothing; percentile clip; 2-means threshold.
    vals = [];
    for i = 1:numel(traces)
        v = traces{i};
        if isempty(v), continue; end
        v = double(v(:)');

        if doSmooth && smoothWinSec > 0
            w = max(1, round(smoothWinSec / dt));
            v = movmean(v, w, 'Endpoints','shrink');
        end

        vals = [vals, v]; %#ok<AGROW>
        if numel(vals) > maxPts*2
            vals = vals(randperm(numel(vals), maxPts*2));
        end
    end

    if numel(vals) < 200
        thr = NaN; return;
    end
    if numel(vals) > maxPts
        vals = vals(randperm(numel(vals), maxPts));
    end

    lo = prctile(vals, clipPrct(1));
    hi = prctile(vals, clipPrct(2));
    vals = vals(vals >= lo & vals <= hi);

    if numel(vals) < 200
        thr = NaN; return;
    end

    try
        [~, C] = kmeans(vals(:), 2, 'Replicates', 5, 'MaxIter', 200);
        C = sort(C(:));
        thr = mean(C); % midpoint
    catch
        thr = NaN;
    end
end

function hasLick = computeHasLickInWindow(traces, dt, stimOnsetSec, winSec, doSmooth, smoothWinSec, thr, minDurSec)
% Trial has lick if (smoothed) signal > thr for any sample in window (or >=minDurSec contiguous).
    n = numel(traces);
    hasLick = false(n,1);
    for i = 1:n
        v = traces{i};
        if isempty(v), continue; end
        v = double(v(:)');

        if doSmooth && smoothWinSec > 0
            w = max(1, round(smoothWinSec / dt));
            v = movmean(v, w, 'Endpoints','shrink');
        end

        tAbs = (0:(numel(v)-1)) * dt;
        tRel = tAbs - stimOnsetSec;

        idx = (tRel >= winSec(1)) & (tRel <= winSec(2));
        if ~any(idx), continue; end

        b = (v(idx) > thr);

        if minDurSec <= 0
            hasLick(i) = any(b);
        else
            minSamp = max(1, round(minDurSec / dt));
            hasLick(i) = any(conv(double(b), ones(1,minSamp), 'same') >= minSamp);
        end
    end
end

function sel = sample_neurons_across_sessions(sessCellCounts, nTarget)
% Return [nTarget x 2] array of (sessionIndex, cellIndex), sampling sessions proportional to cell counts.
    nSess = numel(sessCellCounts);
    sel = zeros(0,2);
    p = sessCellCounts / sum(sessCellCounts);
    while size(sel,1) < nTarget
        s = randsample(nSess, 1, true, p);
        c = randi(sessCellCounts(s));
        pair = [s, c];
        if ~ismember(pair, sel, 'rows')
            sel(end+1,:) = pair; %#ok<AGROW>
        end
    end
end

function [sessList, colsBySess, nNeurBySess, nNeurTotal, wBySess] = build_session_column_map(sel, nSess)
% For selected neurons sel(:,1)=session, sel(:,2)=cellIndex, build mapping to columns 1..nNeurTotal.
    nNeurBySess = zeros(nSess,1);
    for s = 1:nSess
        nNeurBySess(s) = sum(sel(:,1)==s);
    end
    nNeurTotal = sum(nNeurBySess);
    colsBySess = cell(nSess,1);

    col0 = 0;
    sessList = find(nNeurBySess > 0)';
    for s = sessList
        ns = nNeurBySess(s);
        colsBySess{s} = (col0+1):(col0+ns);
        col0 = col0 + ns;
    end
    wBySess = nNeurBySess; % weights for lick pseudo-trace
end

function Xi = build_Xi_blocks(sessionsT, sel, winIdx)
% Return Xi{session} = [nTrials_s x nNeur_s] window-summed responses for selected neurons in that session.
    nSess = numel(sessionsT);
    Xi = cell(nSess,1);

    for s = 1:nSess
        Ts = sessionsT{s};
        mask = (sel(:,1) == s);
        if ~any(mask)
            Xi{s} = [];
            continue;
        end
        cellsLocal = sel(mask,2);
        Xi{s} = squeeze(nansum(Ts(:, cellsLocal, winIdx), 3)); % [nTrials_s x nSelCells_s]
    end
end

function draw = draw_trial_ids_once(sessionsV, sessionsHasLick, sessList, goVal, nogoVal, mode, nPseudo, minTrialsPerCond, strictPerSessionEligibility)
% Sample trial indices ONCE per replicate+mode; reused across all timepoints.
% trialsGoUsed(p,ii) is a REAL trial index within session s = sessIdxVec(ii).
    sessIdxVec = sessList(:)';
    nSessContrib = numel(sessIdxVec);

    trialsGoUsed   = NaN(nPseudo, nSessContrib);
    trialsNoGoUsed = NaN(nPseudo, nSessContrib);

    poolsGo   = cell(max(sessIdxVec),1);
    poolsNoGo = cell(max(sessIdxVec),1);

    okSess = true(1,nSessContrib);

    for ii = 1:nSessContrib
        s = sessIdxVec(ii);
        V = double(sessionsV{s}(:));
        L = logical(sessionsHasLick{s}(:));

        idxGo   = (V == goVal);
        idxNoGo = (V == nogoVal);

        switch lower(mode)
            case 'all'
                mGo   = idxGo;
                mNoGo = idxNoGo;
            case 'lick'
                mGo   = idxGo   & L;
                mNoGo = idxNoGo & L;
            case 'nolick'
                mGo   = idxGo   & ~L;
                mNoGo = idxNoGo & ~L;
            otherwise
                error('Unknown mode: %s', mode);
        end

        poolsGo{s}   = find(mGo);
        poolsNoGo{s} = find(mNoGo);

        if strictPerSessionEligibility
            if numel(poolsGo{s}) < minTrialsPerCond || numel(poolsNoGo{s}) < minTrialsPerCond
                okSess(ii) = false;
            end
        end
    end

    if strictPerSessionEligibility && ~all(okSess)
        % If strict: drop those sessions entirely by leaving NaNs in draw (their neurons become unusable)
        % but keep structure.
    end

    for p = 1:nPseudo
        for ii = 1:nSessContrib
            if strictPerSessionEligibility && ~okSess(ii)
                continue;
            end
            s = sessIdxVec(ii);

            pg = poolsGo{s};
            pn = poolsNoGo{s};

            if ~isempty(pg), trialsGoUsed(p,ii)   = pg(randi(numel(pg))); end
            if ~isempty(pn), trialsNoGoUsed(p,ii) = pn(randi(numel(pn))); end
        end
    end

    draw = struct('sessIdxVec', sessIdxVec, ...
                  'trialsGoUsed', trialsGoUsed, ...
                  'trialsNoGoUsed', trialsNoGoUsed, ...
                  'mode', mode);
end

function [sessionsLickGrid, tLickGrid] = build_sessionsLickGrid(grp, stimOnsetSec, doSmooth, smoothWinSec)
% Build sessionsLickGrid{s} = [nTrials_s x L] where L is max lick length within that session.
% If lick traces are missing (Response shortcut), this returns empty grids (but aux still works).
    nSess = numel(grp.sessionsT);
    sessionsLickGrid = cell(nSess,1);
    tLickGrid = [];

    for s = 1:nSess
        rawTraces = grp.sessionsLickRaw{s};
        dtByTrial = grp.sessionsLickDt{s};

        if isempty(rawTraces) || all(cellfun(@isempty, rawTraces))
            sessionsLickGrid{s} = [];
            continue;
        end

        dtSess = median(dtByTrial(isfinite(dtByTrial)));
        if ~(isfinite(dtSess) && dtSess > 0)
            sessionsLickGrid{s} = [];
            continue;
        end

        Lmax = max(cellfun(@numel, rawTraces));
        M = NaN(numel(rawTraces), Lmax);

        for i = 1:numel(rawTraces)
            v = rawTraces{i};
            if isempty(v), continue; end
            v = double(v(:)');

            if doSmooth && smoothWinSec > 0
                w = max(1, round(smoothWinSec / dtSess));
                v = movmean(v, w, 'Endpoints','shrink');
            end

            M(i,1:numel(v)) = v;
        end

        sessionsLickGrid{s} = M;

        if isempty(tLickGrid)
            tAbs = (0:(Lmax-1))*dtSess;
            tLickGrid = tAbs - stimOnsetSec;
        end
    end
end

function [lickPseudoGo, lickPseudoNoGo] = build_weighted_pseudo_lick(sessionsLickGrid, draw, wBySess, tLickGrid)
% Weighted average lick trace of each pseudo-trial across contributing sessions.
% Weight = #neurons contributed by that session (wBySess(s)).
% If lick grid missing for some session, that session contributes 0 to lick.
    sessIdxVec = draw.sessIdxVec;
    nPseudo = size(draw.trialsGoUsed,1);

    if isempty(tLickGrid)
        lickPseudoGo = [];
        lickPseudoNoGo = [];
        return;
    end

    L = numel(tLickGrid);
    lickPseudoGo   = zeros(nPseudo, L);
    lickPseudoNoGo = zeros(nPseudo, L);

    Wtot = sum(wBySess(sessIdxVec));
    if Wtot <= 0, Wtot = 1; end

    for p = 1:nPseudo
        accG = zeros(1,L);
        accN = zeros(1,L);

        for ii = 1:numel(sessIdxVec)
            s = sessIdxVec(ii);
            w = wBySess(s);

            M = sessionsLickGrid{s};
            if isempty(M), continue; end

            tiG = draw.trialsGoUsed(p,ii);
            tiN = draw.trialsNoGoUsed(p,ii);

            if isfinite(tiG) && tiG >= 1 && tiG <= size(M,1)
                v = M(tiG,:);
                v(~isfinite(v)) = 0;
                if numel(v) < L, v(1,L) = 0; end %#ok<AGROW>
                accG = accG + w * v(1:L);
            end

            if isfinite(tiN) && tiN >= 1 && tiN <= size(M,1)
                v = M(tiN,:);
                v(~isfinite(v)) = 0;
                if numel(v) < L, v(1,L) = 0; end %#ok<AGROW>
                accN = accN + w * v(1:L);
            end
        end

        lickPseudoGo(p,:)   = accG / Wtot;
        lickPseudoNoGo(p,:) = accN / Wtot;
    end
end

function cap = capacity_from_fixed_draw(Xi, colsBySess, sessList, goVal, nogoVal, draw, saveFolder, tag, capacitySuffix)
% Build Pgo/Pnogo (raw pseudo-trials) from fixed trial IDs across sessions, then run check_capacity.
    cap = NaN;

    nPseudo = size(draw.trialsGoUsed,1);
    nNeurTotal = 0;
    for s = sessList
        nNeurTotal = max(nNeurTotal, max(colsBySess{s}));
    end
    if nNeurTotal == 0, return; end

    Pgo   = NaN(nPseudo, nNeurTotal);
    Pnogo = NaN(nPseudo, nNeurTotal);

    sessIdxVec = draw.sessIdxVec;

    for p = 1:nPseudo
        for ii = 1:numel(sessIdxVec)
            s = sessIdxVec(ii);
            cols = colsBySess{s};
            if isempty(cols), continue; end

            Xs = Xi{s};
            if isempty(Xs), continue; end

            tiG = draw.trialsGoUsed(p,ii);
            tiN = draw.trialsNoGoUsed(p,ii);

            if isfinite(tiG) && tiG >= 1 && tiG <= size(Xs,1)
                Pgo(p, cols) = Xs(tiG, :);
            end
            if isfinite(tiN) && tiN >= 1 && tiN <= size(Xs,1)
                Pnogo(p, cols) = Xs(tiN, :);
            end
        end
    end

    % If too many NaNs, check_capacity may fail or return NaNs.
    tuning_function = cat(3, Pgo', Pnogo'); %#ok<NASGU>
    data_titles     = {tag};                %#ok<NASGU>
    session_ids     = 1;                    %#ok<NASGU>

    if ~exist(saveFolder,'dir'), mkdir(saveFolder); end
    name_prefix = fullfile(saveFolder, ['super_' tag]);
    save([name_prefix '_tuning.mat'], 'tuning_function', 'data_titles', '-v7.3');

    try
        check_capacity(name_prefix, session_ids);
    catch
        cap = NaN;
        return;
    end

    filePath = [name_prefix capacitySuffix];
    if ~isfile(filePath)
        d = dir(fullfile(saveFolder, ['*' capacitySuffix]));
        if isempty(d), cap = NaN; return;
        else, [~,ix]=max([d.datenum]); filePath = fullfile(d(ix).folder,d(ix).name);
        end
    end

    try
        mf = matfile(filePath);
        cap = 2 ./ mf.capacity_results;
        if numel(cap) > 1, cap = cap(1); end
    catch
        cap = NaN;
    end
end

function plot_group_mean_sem(t, M, col, labelStr)
% Plot mean ± SEM for an [nRep x nTime] matrix
    if isempty(M), return; end
    mu  = nanmean(M,1);
    sem = nanstd(M,0,1) ./ sqrt(sum(~isnan(M),1));
    face = col + (1-col)*0.35;
    face = min(max(face,0),1);
    fill([t fliplr(t)], [mu+sem fliplr(mu-sem)], face, 'EdgeColor','none', 'FaceAlpha',0.25, 'HandleVisibility','off');
    plot(t, mu, '-', 'Color', col, 'LineWidth', 2.6, 'DisplayName', labelStr);
end
