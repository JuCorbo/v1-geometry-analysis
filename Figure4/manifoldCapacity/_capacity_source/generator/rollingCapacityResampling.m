%% Super-mouse capacity timecourse (all days, Naive vs Trained)
% - Time window: -1s to +4s
% - Rolling & expanding windows, step=3 frames, win=3 frames
% - 5 pseudo-sessions per group/day, 500 neurons, 150 pseudotrials/cond
% - Capacity readout: 2 ./ capacity_results from *_capacity_balanced_s1.mat

clear; clc; close all;

%% -------------------- USER SETTINGS --------------------
rootDir    = 'D:/embedding/TraceSpikes/';

dbNaive    = 'BatuhanADNaive';
dbTrained  = 'BatuhanAngDist';      % <-- change if your trained DB string differs

nPseudoSessions       = 5;          % per group per day
nNeuronsTarget        = 500;        % neurons per pseudo-session
nPseudotrialsPerCond  = 150;        % per class (Go, NoGo)

% Time axis / windows
dt        = 0.0646;                 % s per frame
rasterLen = 140;
xval      = linspace(0,(rasterLen-1)*dt,rasterLen) - 1.5;  % 0 = stim onset
tcWindowSec = [-1 +4];              % analysis window
winFrames   = 3;                    % rolling window length (~190 ms)
stepFrames  = 3;                    % stride (~190 ms)

% Baseline subtraction
baselineWindowSec = [-1.5 0];
subtractBaseline  = true;

% Capacity file suffix produced by your toolbox
capacitySuffix = '_capacity_balanced_s1.mat';

% Plot colors
colNaiveRoll   = [0.20 0.40 0.90];
colNaiveExp    = [0.20 0.60 0.90];
colTrainedRoll = [0.90 0.20 0.20];
colTrainedExp  = [0.90 0.40 0.20];

%% -------------------- DERIVED INDICES --------------------
baselineWindow = find(xval >= baselineWindowSec(1) & xval <= baselineWindowSec(2));
tcIdx     = find(xval >= tcWindowSec(1) & xval <= tcWindowSec(2));
startIdx  = tcIdx(1);
endIdx    = tcIdx(end);
centers   = startIdx:stepFrames:endIdx;
tvec      = xval(centers);

%% -------------------- DAY LOOP --------------------
dayFolders = dir(fullfile(rootDir, 'D*'));
dayFolders = dayFolders([dayFolders.isdir]);

if isempty(dayFolders)
    error('No day folders (D*) found under %s', rootDir);
end

for d = 1:numel(dayFolders)
    day = dayFolders(d).name;
    fprintf('\n\n========================= %s =========================\n', day);

    outRoot = fullfile(rootDir, day, 'supermouse');
    if ~exist(outRoot,'dir'), mkdir(outRoot); end

    % Pool traces across all sessions for this day & group
    [groupNaive,  nCellsNaive ] = load_group_pool(rootDir, day, dbNaive,  rasterLen, subtractBaseline, baselineWindow, endIdx);
    [groupTrained,nCellsTrained] = load_group_pool(rootDir, day, dbTrained,rasterLen, subtractBaseline, baselineWindow, endIdx);

    % Compute pseudo-sessions (super-mouse) and capacity timecourses
    fprintf('\n=== NAIVE (%s) ===\n', dbNaive);
    resNaive   = run_supermouse_group(groupNaive,   nPseudoSessions, nNeuronsTarget, ...
                                      centers, tvec, winFrames, startIdx, capacitySuffix, outRoot, sprintf('%s_Naive', day), nPseudotrialsPerCond);

    fprintf('\n=== TRAINED (%s) ===\n', dbTrained);
    resTrained = run_supermouse_group(groupTrained, nPseudoSessions, nNeuronsTarget, ...
                                      centers, tvec, winFrames, startIdx, capacitySuffix, outRoot, sprintf('%s_Trained', day), nPseudotrialsPerCond);

    % Plot summary per day
    fig = figure('Color','w','Position',[120 120 980 560]); hold on;
    plot_group_mean_sem(tvec, resNaive.capacityRolling,   colNaiveRoll,   'Naive (rolling)');
    plot_group_mean_sem(tvec, resNaive.capacityExpanding, colNaiveExp,    'Naive (expanding)', true);
    plot_group_mean_sem(tvec, resTrained.capacityRolling, colTrainedRoll, 'Trained (rolling)');
    plot_group_mean_sem(tvec, resTrained.capacityExpanding,colTrainedExp, 'Trained (expanding)', true);
    xline(0,'--k','Stim onset','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','left');
    xlabel('Time from stimulus onset (s)');
    ylabel('Capacity (2 / capacity\_results)');
    title(sprintf('Super-mouse — %s — %s vs %s | %d neurons, %d pseudotrials/cond, %d pseudo-sessions', ...
          day, dbNaive, dbTrained, nNeuronsTarget, nPseudotrialsPerCond, nPseudoSessions), 'Interpreter','none');
    legend('Location','best'); grid on; box off;

    saveas(fig, fullfile(outRoot, sprintf('supermouse_%s_summary.png', day)));
    save(fullfile(outRoot, sprintf('supermouse_%s_summary.mat', day)), ...
         'resNaive','resTrained','tvec','day','dbNaive','dbTrained', ...
         'nNeuronsTarget','nPseudotrialsPerCond','nPseudoSessions','winFrames','stepFrames');
    close(fig);

    fprintf('Saved outputs to %s\n', outRoot);
end

disp('All days processed.');

%% ================== FUNCTIONS ==================

function [grp, totalCells] = load_group_pool(rootDir, day, dbName, rasterLen, subtractBaseline, baselineWindow, endIdx)
% Gather per-session traces for a given day & DB by re-querying SQL.
% Reconstruct T(trial,cell,time) and Vstim(trial,1) per session.
    dayDir = fullfile(rootDir, day, 'timecourse');
    if ~exist(dayDir,'dir'); warning('Day folder not found: %s', dayDir); grp = struct('sessionsT',{{}},'sessionsV',{{}}); totalCells=0; return; end

    summaries = dir(fullfile(dayDir, sprintf('*_%s_capacity_timecourse.mat', dbName)));
    if isempty(summaries)
        warning('No session summaries for %s %s. Pool will be empty.', day, dbName);
        grp = struct('sessionsT',{{}},'sessionsV',{{}});
        totalCells = 0;
        return;
    end

    conn = DBconnection({dbName});

    sessionsT = {}; sessionsV = {}; totalCells = 0;

    for k = 1:numel(summaries)
        sessName = extractBefore(summaries(k).name, ['_' dbName '_capacity_timecourse.mat']); % session id
        selectquery = [' select ctt.Experiment,' ...
            ' ctt.Cell ,ctt.Trial, ctt.TraceSpikes, tt.Visual_Stim ' ...
            ' from Cell_Trial_Table AS ctt ' ...
            ' inner join Trial_Table tt ' ...
            ' on tt.Experiment = ctt.Experiment AND tt.Trial=ctt.Trial' ...
            ' inner join Cell_Table ct ' ...
            ' on ct.Experiment = ctt.Experiment AND ct.Cell = ctt.Cell ' ...
            ' WHERE ct.isCell > 0.8 ' ...
            ' AND ct.DeconvCorr > 0.8 ' ...
            ' AND ct.spikeProb_TCTrials > 0.1' ...
            ' AND ct.roundness > 0.2' ...
            ' AND ct.nPix*ct.magFactor > 15' ...
            ' AND (tt.Block != ''Orientation Tuning'') ' ...
            ' AND ct.Experiment= ''' sessName '''' ];

        S = select(conn, selectquery);
        if isempty(S), warning('No task data for %s (%s)', sessName, dbName); continue; end
        S = sortrows(S, {'Trial','Cell'});

        unique_trials = unique(S.Trial);
        unique_cells  = unique(S.Cell);
        nTrials = numel(unique_trials);
        nCells  = numel(unique_cells);

        T = NaN(nTrials, nCells, rasterLen, 'double');
        Vstim = NaN(nTrials,1);

        for i = 1:nTrials
            rows = S.Trial == unique_trials(i);
            trial_cells = S.Cell(rows);
            rawTraces   = S.TraceSpikes(rows);
            Vstim(i)    = unique(S.Visual_Stim(rows));
            [~, loc]    = ismember(trial_cells, unique_cells);
            for j = 1:numel(loc)
                tr = getArrayFromByteStream(typecast(rawTraces{j}, 'uint8'));
                tr = double(tr(:)');
                tr = tr(1:min(end, rasterLen));
                if numel(tr) < endIdx, continue; end
                if subtractBaseline
                    tr = tr - mean(tr(baselineWindow),'omitnan');
                end
                T(i, loc(j), 1:numel(tr)) = tr;
            end
        end

        sessionsT{end+1} = T; %#ok<AGROW>
        sessionsV{end+1} = Vstim; %#ok<AGROW>
        totalCells = totalCells + nCells;
    end

    grp = struct('sessionsT', {sessionsT}, 'sessionsV', {sessionsV});
end

function res = run_supermouse_group(grp, nRep, nNeuronsTarget, centers, tvec, winFrames, startIdx, capacitySuffix, outRoot, tagBase, nPseudotrialsPerCond)
% Build pseudo-sessions and compute capacity timecourses; returns matrices [nRep x nTime]
    if isempty(grp.sessionsT)
        warning('Group %s has no sessions. Returning NaNs.', tagBase);
        res = struct('capacityRolling', NaN(nRep,numel(tvec)), ...
                     'capacityExpanding', NaN(nRep,numel(tvec)), ...
                     'tvec', tvec, 'group', tagBase);
        return;
    end

    saveFolder = fullfile(outRoot, ['super_' tagBase]);
    if ~exist(saveFolder,'dir'), mkdir(saveFolder); end

    % Build time windows
    rollingWindows  = arrayfun(@(t) max(startIdx, t-winFrames+1):t, centers, 'uni', 0);
    expandingWindow = arrayfun(@(t) startIdx:t, centers, 'uni', 0);

    capacityRolling   = NaN(nRep, numel(centers));
    capacityExpanding = NaN(nRep, numel(centers));

    % Session-wise cell counts for neuron sampling
    sessCellCounts = cellfun(@(Ts) size(Ts,2), grp.sessionsT);

    for r = 1:nRep
        fprintf('  Pseudo-session %d/%d...\n', r, nRep);

        % Sample neurons across sessions to reach target count
        sel = sample_neurons_across_sessions(sessCellCounts, nNeuronsTarget);

        for k = 1:numel(centers)
            % Rolling
            [Xr, Vr] = sum_window_supermouse(grp.sessionsT, grp.sessionsV, sel, rollingWindows{k});
            capR = capacity_from_matrix(Xr, Vr, saveFolder, sprintf('%s_roll_t%03d_rep%02d', tagBase, centers(k), r), capacitySuffix, nPseudotrialsPerCond);
            capacityRolling(r, k) = capR;

            % Expanding
            [Xe, Ve] = sum_window_supermouse(grp.sessionsT, grp.sessionsV, sel, expandingWindow{k});
            capE = capacity_from_matrix(Xe, Ve, saveFolder, sprintf('%s_expand_t%03d_rep%02d', tagBase, centers(k), r), capacitySuffix, nPseudotrialsPerCond);
            capacityExpanding(r, k) = capE;
        end
    end

    res = struct('capacityRolling', capacityRolling, ...
                 'capacityExpanding', capacityExpanding, ...
                 'tvec', tvec, 'group', tagBase);
end

function sel = sample_neurons_across_sessions(sessCellCounts, nTarget)
% Return [nTarget x 2] array of (sessionIndex, cellIndex), sampling sessions
% proportional to their cell counts, with uniqueness enforced.
    nSess = numel(sessCellCounts);
    sel = zeros(0, 2);
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

function [X, V] = sum_window_supermouse(sessionsT, sessionsV, selCells, winIdx)
% Build trial-by-neuron response matrix across all sessions using selected neurons.
% Output:
%   X: [nTrialsAll x nNeurons]
%   V: [nTrialsAll x 1]
    Xblocks = cell(numel(sessionsT),1);
    Vblocks = cell(numel(sessionsV),1);

    for s = 1:numel(sessionsT)
        Ts = sessionsT{s};
        Vs = sessionsV{s};
        mask = (selCells(:,1) == s);
        if ~any(mask)
            Xblocks{s} = []; Vblocks{s} = [];
            continue;
        end
        cellsLocal = selCells(mask,2);
        Xi = squeeze(nansum(Ts(:, cellsLocal, winIdx), 3)); % [nTrials_s x nSelCells_s]
        Xblocks{s} = Xi;
        Vblocks{s} = Vs;
    end

    X = cat_rows_with_fill(Xblocks);
    V = cat_rows_with_fill(Vblocks);
end

function out = cat_rows_with_fill(blocks)
% Concatenate vertically; pad columns to the max if needed.
    if isempty(blocks), out = []; return; end
    nonempty = ~cellfun(@isempty, blocks);
    blocks = blocks(nonempty);
    if isempty(blocks), out = []; return; end
    nCols = max(cellfun(@(x) size(x,2), blocks));
    for i = 1:numel(blocks)
        if size(blocks{i},2) < nCols
            blocks{i}(:, end+1:nCols) = NaN;
        end
    end
    out = vertcat(blocks{:});
end

function cap = capacity_from_matrix(X, Vstim, saveFolder, tag, capacitySuffix, nPseudotrialsPerCond)
% Create pseudotrials per condition, run check_capacity, read 2./capacity_results.
    cap = NaN;

    isGo   = (Vstim == 45);
    isNoGo = ~isGo;
    Xgo    = X(isGo, :);
    Xnogo  = X(isNoGo, :);

    if isempty(Xgo) || isempty(Xnogo)
        warning('Insufficient trials per condition for %s', tag);
        return;
    end

    % Group sizes for bootstrap-averaging to reach exactly nPseudotrialsPerCond
    gGo   = max(1, ceil(size(Xgo,1)   / nPseudotrialsPerCond));
    gNoGo = max(1, ceil(size(Xnogo,1) / nPseudotrialsPerCond));

    Pgo   = make_pseudotrials(Xgo,   nPseudotrialsPerCond, gGo);   % [nP x nNeurons]
    Pnogo = make_pseudotrials(Xnogo, nPseudotrialsPerCond, gNoGo);

    % Assemble tuning_function [nNeurons x nTrials x 2]
    tuning_function = cat(3, Pgo', Pnogo'); %#ok<NASGU>
    data_titles     = {tag};                %#ok<NASGU>
    session_ids     = 1;                    %#ok<NASGU>

    name_prefix = fullfile(saveFolder, ['super_' tag]);
    save([name_prefix '_tuning.mat'], 'tuning_function', 'data_titles', '-v7.3');

    try
        check_capacity(name_prefix, session_ids);
    catch ME
        warning('check_capacity failed (%s): %s', tag, ME.message);
        return;
    end

    filePath = [name_prefix capacitySuffix];
    if ~isfile(filePath)
        % Fallback: try newest matching capacity file in folder
        d = dir(fullfile(saveFolder, ['*' capacitySuffix]));
        if isempty(d), return; else, [~,ix]=max([d.datenum]); filePath = fullfile(d(ix).folder,d(ix).name); end
    end

    try
        mf = matfile(filePath);
        cap = 2 ./ mf.capacity_results;   % expected scalar
        if numel(cap) > 1, cap = cap(1); end
    catch ME
        warning('Failed reading capacity_results in %s: %s', filePath, ME.message);
        cap = NaN;
    end
end

function P = make_pseudotrials(X, nPseudo, groupSize)
% Create nPseudo pseudotrials by averaging "groupSize" randomly sampled trials (with replacement).
% X: [nTrials x nNeurons] -> P: [nPseudo x nNeurons]
    nTrials = size(X,1);
    P = NaN(nPseudo, size(X,2));
    for i = 1:nPseudo
        idx = randi(nTrials, groupSize, 1);
        P(i,:) = mean(X(idx, :), 1, 'omitnan');
    end
end

function plot_group_mean_sem(t, M, col, labelStr, altShade)
% Plot mean ± SEM for an [nRep x nTime] matrix
    if nargin < 5, altShade = false; end
    if isempty(M), return; end
    mu  = nanmean(M,1);
    sem = nanstd(M,0,1) ./ sqrt(sum(~isnan(M),1));
    face = col + (1-col)*(0.50*(altShade) + 0.30*(~altShade));
    face = min(max(face,0),1);
    fill([t fliplr(t)], [mu+sem fliplr(mu-sem)], face, 'EdgeColor','none', 'FaceAlpha',0.25, 'DisplayName',[labelStr ' ± SEM']);
    plot(t, mu, '-', 'Color', col, 'LineWidth', 2.4, 'DisplayName', labelStr);
end
