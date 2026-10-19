%% ============================================================
%  rollingCapacity_randomThirdsAvg.m
%
%  Goal:
%   For each session in BOTH databases (naive + trained),
%   randomly split trials into 3 equal parts (stratified within Go/NoGo),
%   compute ROLLING capacity timecourse for each split,
%   then average the 3 timecourses -> one "third-matched" timecourse per session.
%
%  Notes:
%   - Capacity logic and QC filters are kept identical to rollingCapacity.m
%   - We query ALL task trials once per session, decode TraceSpikes once,
%     then index into T/Vstim for each split (no repeated SQL queries).
%   - Stratified split avoids missing Go or NoGo in a split.
%
%  Output:
%   D:/embedding/TraceSpikes/D{day}/timecourse/
%     <session>_<DB>_randThirdAvg_capacity_timecourse.mat
%   also saves the 3 split timecourses (for sanity checks).
%% ============================================================

clear; clc; close all;

%% ---------------------- USER PARAMS ----------------------
DBs = {'BatuhanADNaive','BatuhanAngDist'};   % naive + trained
BehavCond = {'D1','D2','D3','D4','D5','D6'};

% TraceSpikes / time axis
dt        = 0.0646;                 % seconds per frame
rasterLen = 140;                    % frames per trial
xval      = linspace(0,(rasterLen-1)*dt,rasterLen) - 1.5;  % 0 = stim onset

% Baseline subtraction
baselineWindowSec = [-1.5 0];
subtractBaseline  = true;

% Timecourse window and step
tcWindowSec = [-1 4];
winFrames   = 3;                    % ~190 ms window length
stepFrames  = 3;                    % ~190 ms stride

% Output root
saveRoot = 'D:/embedding/TraceSpikes/';

% Capacity filename pattern produced by check_capacity
capacity_suffix = '_capacity_balanced_s1.mat';

% Randomness
rngSeedBase = 1;  % deterministic across runs; per-session seed is derived from this

%% ---------------------- DERIVED INDICES ----------------------
baselineWindow = find(xval >= baselineWindowSec(1) & xval <= baselineWindowSec(2));
tcIdx     = find(xval >= tcWindowSec(1) & xval <= tcWindowSec(2));
startIdx  = tcIdx(1);
endIdx    = tcIdx(end);
centers   = startIdx:stepFrames:endIdx;
tvec      = xval(centers);

rollingWindows = arrayfun(@(t) max(startIdx, t-winFrames+1):t, centers, 'uni', 0);

%% ---------------------- MAIN LOOPS ----------------------
for d = 1:numel(DBs)
    DB = DBs{d};
    fprintf('\n============================================================\n');
    fprintf('DB: %s\n', DB);

    conn = DBconnection({DB});

    % Pref-orientation field depends on DB (same logic as your scripts)
    if strcmp(DB, 'BatuhanAngDist')
        stringPrefOri = 'Pref_Orientation_spikes_2';
    elseif strcmp(DB, 'BatuhanADNaive')
        stringPrefOri = 'Pref_Orientation_spikes_1';
    else
        error('Unknown DB name: %s', DB);
    end

    for exps = 1:numel(BehavCond)
        tic;
        fprintf('\n------------------------------------------------------------\n');
        fprintf('Processing day %d (%s)\n', exps, BehavCond{exps});

        % Session list (UNCHANGED)
        selectquery = [' select distinct ct.Experiment' ...
            ' from Cell_Table AS ct ' ...
            ' inner join Trial_Table tt ' ...
            ' on tt.Experiment = ct.Experiment ' ...
            ' WHERE ct.' stringPrefOri ' IS NOT NULL ' ...
            ' AND ct.isCell > 0.8 ' ...
            ' AND ct.DeconvCorr > 0.8 ' ...
            ' AND ct.spikeProb_TCTrials > 0.1' ...
            ' AND ct.roundness > 0.2' ...
            ' AND ct.nPix*ct.magFactor > 15' ...
            ' AND tt.Behav_Cond = ''' BehavCond{exps} '''' ];

        allExps = select(conn, selectquery);
        if isempty(allExps)
            warning('No sessions for %s %s', DB, BehavCond{exps});
            continue;
        end

        saveFolderDay = fullfile(saveRoot, ['D' num2str(exps)], 'timecourse');
        if ~exist(saveFolderDay, 'dir'), mkdir(saveFolderDay); end

        for sess = 1:height(allExps)
            sessName = allExps.Experiment{sess};
            fprintf('\n=== %s — Day %d/%d — Session %d/%d (%s) ===\n', ...
                DB, exps, numel(BehavCond), sess, height(allExps), sessName);

            % ----------- Query trial-by-cell traces (UNCHANGED) -----------
            selectquery = [' select ctt.Experiment,' ...
                ' ctt.Cell ,ctt.Trial, ctt.TraceSpikes, tt.Visual_Stim ' ...
                ' from Cell_Trial_Table AS ctt ' ...
                ' inner join Trial_Table tt ' ...
                ' on tt.Experiment = ctt.Experiment AND tt.Trial=ctt.Trial' ...
                ' inner join Cell_Table ct ' ...
                ' on ct.Experiment = ctt.Experiment AND ct.Cell = ctt.Cell ' ...
                ' WHERE ct.' stringPrefOri ' IS NOT NULL ' ...
                ' AND (tt.Block != ''Orientation Tuning'' ) ' ...
                ' AND ct.isCell > 0.8 ' ...
                ' AND ct.DeconvCorr > 0.8 ' ...
                ' AND ct.spikeProb_TCTrials > 0.1' ...
                ' AND ct.roundness > 0.2' ...
                ' AND ct.nPix*ct.magFactor > 15' ...
                ' AND tt.Behav_Cond = ''' BehavCond{exps} '''' ...
                ' AND ct.Experiment= ''' sessName '''' ];

            taskData = sortrows(select(conn, selectquery), {'Trial','Cell'});
            if isempty(taskData)
                warning('Empty taskData for %s', sessName);
                continue;
            end

            unique_trials = unique(taskData.Trial);
            unique_cells  = unique(taskData.Cell);
            nTrials = numel(unique_trials);
            nCells  = numel(unique_cells);

            % ----------- Decode & baseline-center traces -----------
            T     = NaN(nTrials, nCells, rasterLen, 'double');
            Vstim = NaN(nTrials, 1);

            for i = 1:nTrials
                rows = taskData.Trial == unique_trials(i);
                trial_cells = taskData.Cell(rows);
                rawTraces   = taskData.TraceSpikes(rows);
                Vstim(i)    = unique(taskData.Visual_Stim(rows));
                [~, loc]    = ismember(trial_cells, unique_cells);

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

            % Need both classes present overall
            isGoAll   = (Vstim == 45);
            isNoGoAll = ~isGoAll;
            if ~any(isGoAll) || ~any(isNoGoAll)
                warning('Skipping %s: missing Go or NoGo trials.', sessName);
                continue;
            end

            % ----------- Random split into thirds (STRATIFIED within Go/NoGo) -----------
            % Use a deterministic per-session seed so reruns are identical
            sessHash = sum(double(char(sessName))) + 1000*exps + 100000*d;
            rng(rngSeedBase + sessHash, 'twister');

            goIdx   = find(isGoAll);
            nogoIdx = find(isNoGoAll);

            goIdx   = goIdx(randperm(numel(goIdx)));
            nogoIdx = nogoIdx(randperm(numel(nogoIdx)));

            goParts   = split_into_thirds(goIdx);
            nogoParts = split_into_thirds(nogoIdx);

            % For each split, combine the corresponding Go+NoGo trials
            cap_roll_splits = NaN(numel(centers), 3);

            for s = 1:3
                idxSplit = [goParts{s}; nogoParts{s}];
                idxSplit = idxSplit(:);

                if ~any(Vstim(idxSplit)==45) || ~any(Vstim(idxSplit)~=45)
                    warning('Split %d in %s lacks Go/NoGo after stratification? Skipping split.', s, sessName);
                    continue;
                end

                Tsub = T(idxSplit, :, :);
                Vsub = Vstim(idxSplit);

                for k = 1:numel(centers)
                    name_roll = fullfile(saveFolderDay, sprintf('%s_%s_randThird_s%d_roll_t%03d', ...
                        sessName, DB, s, centers(k)));

                    cap_roll_splits(k, s) = capacity_from_window(Tsub, Vsub, rollingWindows{k}, name_roll, capacity_suffix);
                end
            end

            % Average across the 3 splits (per timepoint)
            capacity_rolling_randThirdAvg = mean(cap_roll_splits, 2, 'omitnan');

            % ----------- Save summary -----------
            meta = struct( ...
                'DB',DB,'session',sessName,'BehavCond',BehavCond{exps}, ...
                'dt',dt,'rasterLen',rasterLen,'subtractBaseline',subtractBaseline, ...
                'baselineWindowSec',baselineWindowSec,'tcWindowSec',tcWindowSec, ...
                'winFrames',winFrames,'stepFrames',stepFrames,'centers_idx',centers, ...
                'rngSeedBase',rngSeedBase,'rngSeedUsed',rngSeedBase + sessHash, ...
                'splitType','stratified_by_visual_stim','nTrials_total',nTrials, ...
                'nTrials_go_total',sum(isGoAll),'nTrials_nogo_total',sum(isNoGoAll));

            outFile = fullfile(saveFolderDay,[sessName '_' DB '_randThirdAvg_capacity_timecourse.mat']);
            save(outFile, ...
                'capacity_rolling_randThirdAvg', 'cap_roll_splits', 'tvec', 'meta', ...
                '-v7.3');

            fprintf('Saved %s\n', outFile);
        end

        fprintf('Finished %s %s in %.1f s\n', DB, BehavCond{exps}, toc);
    end
end

disp('All done.');

%% ---------------------- LOCAL FUNCTIONS ----------------------
function parts = split_into_thirds(idx)
    % idx: vector of indices
    % returns parts{1..3} approximately equal sizes (difference <= 1)
    n = numel(idx);
    n1 = floor(n/3);
    n2 = floor((n-n1)/2);
    n3 = n - n1 - n2;

    parts = cell(3,1);
    parts{1} = idx(1:n1);
    parts{2} = idx(n1+1 : n1+n2);
    parts{3} = idx(n1+n2+1 : n1+n2+n3);
end

function cap = capacity_from_window(T,Vstim,winIdx,name_prefix,capacity_suffix)
    % Build trial x neuron responses for the given window and compute capacity.
    cap = NaN;
    if isempty(winIdx), return; end

    X = squeeze(nansum(T(:,:,winIdx),3));       % [nTrials x nCells]
    isGo = (Vstim==45); isNoGo = ~isGo;
    nG = sum(isGo); nN = sum(isNoGo);
    if nG==0 || nN==0, return; end

    % Balance classes into tuning_function: [nCells x nTrials x 2]
    if nG >= nN
        TF = cat(3, X(isGo,:)', [X(isNoGo,:); NaN(nG-nN, size(X,2))]');
    else
        TF = cat(3, [X(isGo,:); NaN(nN-nG, size(X,2))]', X(isNoGo,:)');
    end
    tuning_function = TF; %#ok<NASGU>
    data_titles     = {name_prefix}; %#ok<NASGU>
    session_ids     = 1; %#ok<NASGU>

    % Save and run capacity
    save([name_prefix '_tuning.mat'], 'tuning_function', 'data_titles', '-v7.3');
    try
        check_capacity(name_prefix, session_ids);
    catch ME
        warning('check_capacity failed at %s: %s', name_prefix, ME.message);
        return;
    end

    % Read capacity as: 2 / capacity_results from the known file
    filePath = [name_prefix capacity_suffix];
    if ~isfile(filePath)
        d = dir([fileparts(name_prefix) filesep '*' capacity_suffix]);
        if isempty(d), return; else, filePath = fullfile(d(1).folder, d(1).name); end
    end

    try
        mf = matfile(filePath);
        cap = 2 ./ mf.capacity_results;  % scalar
        if numel(cap) > 1, cap = cap(1); end
    catch ME
        warning('Failed reading capacity_results in %s: %s', filePath, ME.message);
        cap = NaN;
    end
end
