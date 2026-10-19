%% Capacity timecourse per session (rolling & expanding windows)
% Control analysis: trained animals, days D4-D5, compute capacity separately for:
%   (1) lick trials only (tt.Response == 1)
%   (2) no-lick trials only (tt.Response == 0)
%
% Uses TraceSpikes -> tuning_function -> check_capacity, then reads:
%   file = matfile(filePath); capacity_val = 2 / file.capacity_results;

%% ---------------------- USER PARAMS ----------------------
DB = 'BatuhanAngDist';
conn = DBconnection({DB});

% Pref-orientation field depends on DB
if strcmp(DB, 'BatuhanAngDist')
    stringPrefOri = 'Pref_Orientation_spikes_2';
elseif strcmp(DB, 'BatuhanADNaive')
    stringPrefOri = 'Pref_Orientation_spikes_1';
else
    error('Unknown DB name');
end

% Trained only, best balance days
BehavCond = {'D4','D5'};

% Lick subsets based on Trial_Table.Response (1 = licked in last second)
lickLabels = {'lickOnly','noLickOnly'};
lickValues = [1, 0];

% TraceSpikes / time axis
dt        = 0.0646;                 % seconds per frame
rasterLen = 140;                    % frames per trial
xval      = linspace(0,(rasterLen-1)*dt,rasterLen) - 1.5;  % 0 = stim onset

% Baseline subtraction
baselineWindowSec = [-1.5 0];
subtractBaseline  = true;

% Timecourse window and step
tcWindowSec = [-1 4];              % from -1s to +4s
winFrames   = 3;                    % ~190 ms window length
stepFrames  = 3;                    % ~190 ms stride

% Output root
saveRoot = 'D:/embedding/TraceSpikes/LickNoLick/';

% Capacity filename pattern (exact filename we will generate)
capacity_suffix = '_capacity_balanced_s1.mat';

useParallel = false;

%% ---------------------- DERIVED INDICES ----------------------
baselineWindow = find(xval >= baselineWindowSec(1) & xval <= baselineWindowSec(2));
tcIdx     = find(xval >= tcWindowSec(1) & xval <= tcWindowSec(2));
startIdx  = tcIdx(1);
endIdx    = tcIdx(end);
centers   = startIdx:stepFrames:endIdx;
tvec      = xval(centers);

%% ---------------------- PER DAY / SESSION LOOPS ----------------------
for exps = 1:numel(BehavCond)
    tic;
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Processing day %d (%s)\n', exps, BehavCond{exps});

    % Session list
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
        warning('No sessions for %s', BehavCond{exps});
        continue;
    end

    for sess = 1:height(allExps)
        sessName = allExps.Experiment{sess};
        fprintf('\n=== Day %s — Session %d/%d (%s) ===\n', ...
            BehavCond{exps}, sess, height(allExps), sessName);

        % ----------- Query trial-by-cell traces (includes Response) -----------
        selectquery = [' select ctt.Experiment,' ...
            ' ctt.Cell ,ctt.Trial, ctt.TraceSpikes, tt.Visual_Stim, tt.Response ' ...
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
        Resp  = NaN(nTrials, 1);

        for i = 1:nTrials
            rows = taskData.Trial == unique_trials(i);
            trial_cells = taskData.Cell(rows);
            rawTraces   = taskData.TraceSpikes(rows);

            % trial-wise metadata (should be constant across cells)
            thisStim = unique(taskData.Visual_Stim(rows));
            thisResp = unique(taskData.Response(rows));
            if numel(thisStim) ~= 1 || numel(thisResp) ~= 1
                warning('Non-unique Visual_Stim or Response within trial %d in %s. Skipping trial.', unique_trials(i), sessName);
                continue;
            end
            Vstim(i) = thisStim;
            Resp(i)  = thisResp;

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

        % ----------- Define windows -----------
        rollingWindows  = arrayfun(@(t) max(startIdx, t-winFrames+1):t, centers, 'uni', 0);
        expandingWindow = arrayfun(@(t) startIdx:t, centers, 'uni', 0);

        % ----------- Compute capacities for lick-only and no-lick-only -----------
        for ls = 1:numel(lickValues)
            thisLickVal   = lickValues(ls);
            thisLickLabel = lickLabels{ls};

            % Filter trials by lick status
            keepTrials = (Resp == thisLickVal) & ~isnan(Vstim);
            Tsub     = T(keepTrials,:,:);
            VstimSub = Vstim(keepTrials);

            if isempty(VstimSub) || size(Tsub,1) < 2
                warning('Skipping %s (%s): too few trials after Response==%d filter.', sessName, thisLickLabel, thisLickVal);
                continue;
            end

            % Need both classes present after filtering
            if ~any(VstimSub==45) || ~any(VstimSub~=45)
                warning('Skipping %s (%s): missing Go or NoGo trials after Response==%d filter.', sessName, thisLickLabel, thisLickVal);
                continue;
            end

            % Save folder per day + lick state
            saveFolderDay = fullfile(saveRoot, BehavCond{exps}, 'timecourse', thisLickLabel);
            if ~exist(saveFolderDay, 'dir'), mkdir(saveFolderDay); end

            capacity_rolling   = NaN(numel(centers),1);
            capacity_expanding = NaN(numel(centers),1);

            if useParallel, gcp('nocreate'); end %#ok<UNRCH>

            for k = 1:numel(centers)
                name_roll = fullfile(saveFolderDay, sprintf('%s_%s_%s_roll_t%03d', sessName, DB, thisLickLabel, centers(k)));
                capacity_rolling(k) = capacity_from_window(Tsub, VstimSub, rollingWindows{k}, name_roll, capacity_suffix);

                name_exp  = fullfile(saveFolderDay, sprintf('%s_%s_%s_expand_t%03d', sessName, DB, thisLickLabel, centers(k)));
                capacity_expanding(k) = capacity_from_window(Tsub, VstimSub, expandingWindow{k}, name_exp, capacity_suffix);
            end

            % ----------- Save summary -----------
            meta = struct('DB',DB,'session',sessName,'BehavCond',BehavCond{exps}, ...
                          'dt',dt,'rasterLen',rasterLen,'subtractBaseline',subtractBaseline, ...
                          'baselineWindowSec',baselineWindowSec,'tcWindowSec',tcWindowSec, ...
                          'winFrames',winFrames,'stepFrames',stepFrames,'centers_idx',centers, ...
                          'lickFilterValue',thisLickVal,'lickFilterLabel',thisLickLabel);

            outFile = fullfile(saveFolderDay, sprintf('%s_%s_%s_capacity_timecourse.mat', sessName, DB, thisLickLabel));
            save(outFile,'capacity_rolling','capacity_expanding','tvec','meta','-v7.3');
            fprintf('Saved %s\n', outFile);
        end
    end

    fprintf('Finished %s in %.1f s\n', BehavCond{exps}, toc);
end
disp('All done.');

%% ---------------------- LOCAL FUNCTIONS ----------------------
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
        % Fallback: try first match in folder
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
