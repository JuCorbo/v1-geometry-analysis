%% ---------------- PERMUTATION TESTS: RADIUS & CENTROID DISTANCE ----------------
% Reads the pull output (baseline-z sigma 0.12); shared groupSeparationTest on
% path; onsets reported at window centre (onsetOffsetSec = dt).
codeDir = fileparts(mfilename('fullpath')); addpath(codeDir);
addpath(fullfile(fileparts(fileparts(fileparts(codeDir))),'_sharedWindowedData','code'));  % canonical groupSeparationTest
load(fullfile(fileparts(codeDir),'cache','originalSpaceBasicGeometrySession.mat'))
dt = 0.0646;
figDir = fullfile(fileparts(codeDir),'panels','main');
if ~exist(figDir, 'dir'), mkdir(figDir); end

geomTimingDir = fullfile(fileparts(codeDir),'stats');
if ~exist(geomTimingDir, 'dir'), mkdir(geomTimingDir); end

% Helper: stack a cell array of [1 x nWin] vectors into [nSess x nWin]
stackSessions = @(C) cell2mat(cellfun(@(v) v(:)', ...
    C(~cellfun(@isempty,C))', 'UniformOutput', false));

measuresToTest = {
    'radiusMix',   radius_bySess,   'left';
    'centDist',    centDist_bySess, 'right';
};

statsGeom = struct();

for mi = 1:size(measuresToTest, 1)
    mName   = measuresToTest{mi,1};
    mData   = measuresToTest{mi,2};
    mTail   = measuresToTest{mi,3};

    statsGeom.(mName) = cell(1, numel(BehavCond));

    fprintf('\n==============================\n');
    fprintf('Testing: %s\n', mName);
    fprintf('==============================\n');

    for exps = 1:numel(BehavCond)
        dataN = stackSessions(mData{1, exps});
        dataT = stackSessions(mData{2, exps});

        if isempty(dataN) || isempty(dataT)
            fprintf('  %s skipped (empty data)\n', BehavCond{exps});
            continue;
        end

        statsGeom.(mName){exps} = groupSeparationTest(dataN, dataT, tvec, ...
            'measureName',  sprintf('%s_%s', mName, BehavCond{exps}), ...
            'tail',         mTail,  ...
            'nPerms',       10000,  ...
            'alpha',        0.05,   ...
            'minSustained', 1,      ...
            'minHedgesG',   0.2,    ...
            'tStart',       0,      ...
            'onsetOffsetSec', dt,   ...
            'correction',   'cluster', ...
            'verbose',      true);
    end
end

%% ---------------- SUMMARY TABLES (per day) ----------------
for mi = 1:size(measuresToTest, 1)
    mName = measuresToTest{mi,1};
    fprintf('\n=== %s - onset times ===\n', mName);
    fprintf('%-6s  %-12s  %-8s  %-8s\n', 'Day', 'Onset (s)', 'nSig', 'nSigSust');
    for exps = 1:numel(BehavCond)
        s = statsGeom.(mName){exps};
        if isempty(s), continue; end
        fprintf('%-6s  %-12s  %-8d  %-8d\n', BehavCond{exps}, ...
            num2str(s.onsetTime, '%.3f'), sum(s.sigMask), sum(s.sigMaskSust));
    end
end

%% ---------------- SAVE TIMING (per day) ----------------
geomTiming.measures   = measuresToTest(:,1)';
geomTiming.days       = BehavCond;
geomTiming.tvec       = tvec;
geomTiming.params.alpha        = 0.05;
geomTiming.params.minHedgesG   = 0.2;
geomTiming.params.minSustained = 1;
geomTiming.params.tStart       = 0;
geomTiming.params.nPerms       = 10000;

for mi = 1:size(measuresToTest, 1)
    mName = measuresToTest{mi,1};
    geomTiming.(mName).onsetTime   = NaN(1, numel(BehavCond));
    geomTiming.(mName).nSig        = zeros(1, numel(BehavCond));
    geomTiming.(mName).nSigSust    = zeros(1, numel(BehavCond));
    geomTiming.(mName).qvals       = cell(1, numel(BehavCond));
    geomTiming.(mName).hedgesG     = cell(1, numel(BehavCond));
    geomTiming.(mName).sigMask     = cell(1, numel(BehavCond));
    geomTiming.(mName).sigMaskSust = cell(1, numel(BehavCond));
    geomTiming.(mName).clusterPmin = NaN(1, numel(BehavCond));

    for exps = 1:numel(BehavCond)
        s = statsGeom.(mName){exps};
        if isempty(s), continue; end
        geomTiming.(mName).onsetTime(exps)   = s.onsetTime;
        geomTiming.(mName).nSig(exps)        = sum(s.sigMask);
        geomTiming.(mName).nSigSust(exps)    = sum(s.sigMaskSust);
        geomTiming.(mName).qvals{exps}       = s.qvals;
        geomTiming.(mName).hedgesG{exps}     = s.hedgesG;
        geomTiming.(mName).sigMask{exps}     = s.sigMask;
        geomTiming.(mName).sigMaskSust{exps} = s.sigMaskSust;
        geomTiming.(mName).clusterPmin(exps) = s.clusterPmin;
    end
end

save(fullfile(geomTimingDir, 'geomTiming.mat'), 'geomTiming');
fprintf('\nGeometry timing saved to: %s\n', fullfile(geomTimingDir, 'geomTiming.mat'));

%% ---------------- POOLED TEST: ALL SESSIONS ACROSS DAYS ----------------
% Pools all sessions across all 6 days per group (~42 naive, ~60 trained).
% Corresponds directly to the grand average in the last summary plot (col 3).
% Non-independence within animal is acknowledged (session-level test).

statsGeomPooled = struct();

fprintf('\n==============================\n');
fprintf('Pooled test (all days)\n');
fprintf('==============================\n');

for mi = 1:size(measuresToTest, 1)
    mName = measuresToTest{mi,1};
    mData = measuresToTest{mi,2};
    mTail = measuresToTest{mi,3};

    allN = [];
    allT = [];
    for exps = 1:numel(BehavCond)
        allN = [allN; stackSessions(mData{1, exps})];  %#ok<AGROW>
        allT = [allT; stackSessions(mData{2, exps})];  %#ok<AGROW>
    end
    nSessNaive   = size(allN, 1);
    nSessTrained = size(allT, 1);

    fprintf('\n--- %s: nNaive=%d  nTrained=%d sessions ---\n', mName, nSessNaive, nSessTrained);

    statsGeomPooled.(mName) = groupSeparationTest(allN, allT, tvec, ...
        'measureName',  sprintf('%s_pooled', mName), ...
        'tail',         mTail,    ...
        'correction',   'cluster', ...
        'nPerms',       10000,    ...
        'alpha',        0.05,     ...
        'minSustained', 1,        ...
        'minHedgesG',   0.2,      ...
        'tStart',       0,        ...
        'onsetOffsetSec', dt,     ...
        'verbose',      true);

    % Append pooled fields to geomTiming
    geomTiming.(mName).onsetTimePooled   = statsGeomPooled.(mName).onsetTime;
    geomTiming.(mName).sigMaskPooled     = statsGeomPooled.(mName).sigMask;
    geomTiming.(mName).sigMaskSustPooled = statsGeomPooled.(mName).sigMaskSust;
    geomTiming.(mName).clusterMaskPooled = statsGeomPooled.(mName).clusterMask;
    geomTiming.(mName).clusterPminPooled = statsGeomPooled.(mName).clusterPmin;
    geomTiming.(mName).hedgesGPooled     = statsGeomPooled.(mName).hedgesG;
    geomTiming.(mName).nSessNaivePooled  = nSessNaive;
    geomTiming.(mName).nSessTrainPooled  = nSessTrained;
end

save(fullfile(geomTimingDir, 'geomTiming.mat'), 'geomTiming');
fprintf('\ngeomTiming resaved with pooled fields.\n');

%% ---------------- SUMMARY TABLE: POOLED ----------------
for mi = 1:size(measuresToTest, 1)
    mName = measuresToTest{mi,1};
    s  = statsGeomPooled.(mName);
    nN = geomTiming.(mName).nSessNaivePooled;
    nT = geomTiming.(mName).nSessTrainPooled;
    fprintf('\n=== %s - POOLED (n=%d naive, %d trained sessions) ===\n', mName, nN, nT);
    fprintf('  Onset:                   %s s\n', num2str(s.onsetTime, '%.3f'));
    fprintf('  Sig (cluster+|G|>=0.2):  %d / %d timepoints\n', sum(s.sigMask),     numel(s.sigMask));
    fprintf('  Sig (sustained):         %d / %d timepoints\n', sum(s.sigMaskSust), numel(s.sigMaskSust));
end

%% ---------------- PLOTTING: DAY MEANS +/- SEM + POOLED +/- SEM ----------------
[train_cent_dayMean, train_cent_daySEM] = day_stats(centDist_bySess(2,:), numel(BehavCond));
[naive_cent_dayMean, naive_cent_daySEM] = day_stats(centDist_bySess(1,:), numel(BehavCond));
[train_rad_dayMean,  train_rad_daySEM]  = day_stats(radius_bySess(2,:),   numel(BehavCond));
[naive_rad_dayMean,  naive_rad_daySEM]  = day_stats(radius_bySess(1,:),   numel(BehavCond));

[train_cent_pooledMean, train_cent_pooledSEM] = pooled_stats(centDist_bySess(2,:));
[naive_cent_pooledMean, naive_cent_pooledSEM] = pooled_stats(centDist_bySess(1,:));
[train_rad_pooledMean,  train_rad_pooledSEM]  = pooled_stats(radius_bySess(2,:));
[naive_rad_pooledMean,  naive_rad_pooledSEM]  = pooled_stats(radius_bySess(1,:));

figure('Color','w','Position',[100 100 2100 900]);
tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
sgtitle(tl, 'MIX geometry summaries: day means \pm SEM and pooled \pm SEM', 'Interpreter','none');

alphaBand    = 0.18;
colTrainBand = blend_with_white(colTrain, 0.65);
colNaiveBand = blend_with_white(colNaive, 0.65);
% D1=darkest (alpha=1.0), D6=faintest (alpha=0.3)
dayAlphas    = linspace(1.0, 0.3, numel(BehavCond));

% ===== Row 1: Centroid distance =====
ax = nexttile(tl,1); hold(ax,'on'); box(ax,'off');
title(ax,'Trained: day mean \pm SEM','Interpreter','none');
xlabel(ax,'Time (s)'); ylabel(ax,'Centroid distance (RMS)');
plot_days_with_alpha(ax, tvec, train_cent_dayMean, train_cent_daySEM, colTrain, colTrainBand, alphaBand, dayAlphas, thickLW);
xlim(ax,[tvec(1) tvec(end)]);

ax = nexttile(tl,2); hold(ax,'on'); box(ax,'off');
title(ax,'Naive: day mean \pm SEM','Interpreter','none');
xlabel(ax,'Time (s)'); ylabel(ax,'Centroid distance (RMS)');
plot_days_with_alpha(ax, tvec, naive_cent_dayMean, naive_cent_daySEM, colNaive, colNaiveBand, alphaBand, dayAlphas, thickLW);
xlim(ax,[tvec(1) tvec(end)]);

% tile 3: difference per day - plotted first so we can use group ylim for scaling
ax_cent_diff = nexttile(tl,3); hold(ax_cent_diff,'on'); box(ax_cent_diff,'off');
title(ax_cent_diff,'Group difference per day (trained \minus naive)','Interpreter','none');
xlabel(ax_cent_diff,'Time (s)'); ylabel(ax_cent_diff,'\Delta Centroid distance (RMS)');
yline(ax_cent_diff, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
for exps = 1:numel(BehavCond)
    mT = train_cent_dayMean(exps,:);
    mN = naive_cent_dayMean(exps,:);
    if all(~isfinite(mT)) || all(~isfinite(mN)), continue; end
    plot(ax_cent_diff, tvec, mT - mN, '-', 'Color', [0 0 0 dayAlphas(exps)], 'LineWidth', thickLW);
end
xlim(ax_cent_diff,[tvec(1) tvec(end)]);
legend(ax_cent_diff, BehavCond, 'Location','best', 'Box','off');

ax = nexttile(tl,4); hold(ax,'on'); box(ax,'off');
title(ax,'Pooled all sessions: trained vs naive \pm SEM','Interpreter','none');
xlabel(ax,'Time (s)'); ylabel(ax,'Centroid distance (RMS)');
plot_mean_with_sem(ax, tvec, train_cent_pooledMean, train_cent_pooledSEM, colTrain, colTrainBand, alphaBand, thickLW);
plot_mean_with_sem(ax, tvec, naive_cent_pooledMean, naive_cent_pooledSEM, colNaive, colNaiveBand, alphaBand, thickLW);
xlim(ax,[tvec(1) tvec(end)]);
legend(ax, {'trained','naive'}, 'Location','best', 'Box','off');
s = statsGeomPooled.centDist;
sigT = s.tvec(s.sigMaskSust);
if ~isempty(sigT)
    yDot = ax.YLim(2) * 0.93;
    scatter(ax, sigT, repmat(yDot, size(sigT)), 6, 'k', 'filled', 'HandleVisibility','off');
end

% Match y-lims for tiles 1,2,4; then set diff plot to same span centered on 0
yl = [inf -inf];
for k = [1 2 4], axk = nexttile(tl,k); yl = [min(yl(1),axk.YLim(1)), max(yl(2),axk.YLim(2))]; end
for k = [1 2 4], ylim(nexttile(tl,k), yl); end
halfSpan = (yl(2) - yl(1)) / 2;
ylim(ax_cent_diff, [-halfSpan halfSpan]);

% ===== Row 2: RadiusMix =====
ax = nexttile(tl,5); hold(ax,'on'); box(ax,'off');
title(ax,'Trained: day mean \pm SEM','Interpreter','none');
xlabel(ax,'Time (s)'); ylabel(ax,'Radius (RMS) [MIX]');
plot_days_with_alpha(ax, tvec, train_rad_dayMean, train_rad_daySEM, colTrain, colTrainBand, alphaBand, dayAlphas, thickLW);
xlim(ax,[tvec(1) tvec(end)]);

ax = nexttile(tl,6); hold(ax,'on'); box(ax,'off');
title(ax,'Naive: day mean \pm SEM','Interpreter','none');
xlabel(ax,'Time (s)'); ylabel(ax,'Radius (RMS) [MIX]');
plot_days_with_alpha(ax, tvec, naive_rad_dayMean, naive_rad_daySEM, colNaive, colNaiveBand, alphaBand, dayAlphas, thickLW);
xlim(ax,[tvec(1) tvec(end)]);

% tile 7: difference per day
ax_rad_diff = nexttile(tl,7); hold(ax_rad_diff,'on'); box(ax_rad_diff,'off');
title(ax_rad_diff,'Group difference per day (trained \minus naive)','Interpreter','none');
xlabel(ax_rad_diff,'Time (s)'); ylabel(ax_rad_diff,'\Delta Radius (RMS) [MIX]');
yline(ax_rad_diff, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
for exps = 1:numel(BehavCond)
    mT = train_rad_dayMean(exps,:);
    mN = naive_rad_dayMean(exps,:);
    if all(~isfinite(mT)) || all(~isfinite(mN)), continue; end
    plot(ax_rad_diff, tvec, mT - mN, '-', 'Color', [0 0 0 dayAlphas(exps)], 'LineWidth', thickLW);
end
xlim(ax_rad_diff,[tvec(1) tvec(end)]);
legend(ax_rad_diff, BehavCond, 'Location','best', 'Box','off');

ax = nexttile(tl,8); hold(ax,'on'); box(ax,'off');
title(ax,'Pooled all sessions: trained vs naive \pm SEM','Interpreter','none');
xlabel(ax,'Time (s)'); ylabel(ax,'Radius (RMS) [MIX]');
plot_mean_with_sem(ax, tvec, train_rad_pooledMean, train_rad_pooledSEM, colTrain, colTrainBand, alphaBand, thickLW);
plot_mean_with_sem(ax, tvec, naive_rad_pooledMean, naive_rad_pooledSEM, colNaive, colNaiveBand, alphaBand, thickLW);
xlim(ax,[tvec(1) tvec(end)]);
legend(ax, {'trained','naive'}, 'Location','best', 'Box','off');
s = statsGeomPooled.radiusMix;
sigT = s.tvec(s.sigMaskSust);
if ~isempty(sigT)
    yDot = ax.YLim(2) * 0.93;
    scatter(ax, sigT, repmat(yDot, size(sigT)), 6, 'k', 'filled', 'HandleVisibility','off');
end

% Match y-lims for tiles 5,6,8; then set diff plot to same span centered on 0
yl = [inf -inf];
for k = [5 6 8], axk = nexttile(tl,k); yl = [min(yl(1),axk.YLim(1)), max(yl(2),axk.YLim(2))]; end
for k = [5 6 8], ylim(nexttile(tl,k), yl); end
halfSpan = (yl(2) - yl(1)) / 2;
ylim(ax_rad_diff, [-halfSpan halfSpan]);

% ===== Vertical reference lines (added after all ylims are final) =====
stimOn  = 0;
stimOff = 3;
% first time bin center (beginning of first window)
firstBin = tvec(1);

colStim  = [0.4 0.4 0.4];  % grey for stim lines
colOnset = [0.8 0.2 0.2];  % red for separation onset

% stim onset/offset on all 8 tiles
for k = 1:8
    axk = nexttile(tl, k);
    xline(axk, stimOn,  '--', 'Color', colStim,  'LineWidth', 0.8, 'HandleVisibility','off');
    xline(axk, stimOff, '--', 'Color', colStim,  'LineWidth', 0.8, 'HandleVisibility','off');
end

% separation onset on pooled panels (tiles 4 and 8)
onCent = statsGeomPooled.centDist.onsetTime;
onRad  = statsGeomPooled.radiusMix.onsetTime;
if ~isnan(onCent)
    xline(nexttile(tl,4), onCent, '-', 'Color', colOnset, 'LineWidth', 1.2, 'HandleVisibility','off');
end
if ~isnan(onRad)
    xline(nexttile(tl,8), onRad,  '-', 'Color', colOnset, 'LineWidth', 1.2, 'HandleVisibility','off');
end
ylim(ax_rad_diff, [-halfSpan halfSpan]);

save_figure_both(gcf, figDir, 'main_centDist_radius_summary');

%% ---------------- PANEL: PEAK centroid distance across days (boxplot + swarm) ----------------
gngAngles = [90 45 30 25 20 15];
postIdx   = tvec >= 0;
peakByDay = cell(2, numel(BehavCond));   % row1 naive, row2 trained
for gi = 1:2
    for d = 1:numel(BehavCond)
        C = centDist_bySess{gi, d};
        C = C(~cellfun(@isempty, C));
        peakByDay{gi, d} = cellfun(@(v) max(v(postIdx)), C);
    end
end
% Kruskal-Wallis: effect of day within each group
allN=[]; gN=[]; allT=[]; gT=[];
for d = 1:numel(BehavCond)
    allN=[allN peakByDay{1,d}]; gN=[gN d*ones(1,numel(peakByDay{1,d}))]; %#ok<AGROW>
    allT=[allT peakByDay{2,d}]; gT=[gT d*ones(1,numel(peakByDay{2,d}))]; %#ok<AGROW>
end
pKW_N = kruskalwallis(allN, gN, 'off');
pKW_T = kruskalwallis(allT, gT, 'off');
fprintf('\n[PEAK centDist day effect] Naive KW p=%.4g | Trained KW p=%.4g\n', pKW_N, pKW_T);
% figure
peakByDayPanel(peakByDay, gngAngles, colNaive, colTrain, ...
    'Peak centroid distance (RMS)', ...
    sprintf('Peak centroid distance by day (KW day effect: Naive p=%.3g, Trained p=%.3g)', pKW_N, pKW_T), ...
    figDir, 'main_peakCentroidDistance_byDay');
% export the per-session peak values + KW to stats/data
try
    fidp = fopen(fullfile(geomTimingDir,'peakCentroidDistance_byDay_stats.txt'),'w');
    fprintf(fidp,'Peak centroid distance (max over t>=0, RMS-norm.), per session, by day.\n');
    fprintf(fidp,'Kruskal-Wallis effect of day within group: Naive p=%.4g ; Trained p=%.4g\n\n', pKW_N, pKW_T);
    fprintf(fidp,'%-8s %-8s %-8s %-8s %-6s\n','group','day','deg','median','n');
    for gi=1:2
        gname = {'Naive','Trained'}; gname=gname{gi};
        for d=1:numel(BehavCond)
            v=peakByDay{gi,d};
            fprintf(fidp,'%-8s %-8s %-8d %-8.3f %-6d\n', gname, BehavCond{d}, gngAngles(d), median(v), numel(v));
        end
    end
    fclose(fidp);
catch
end

%% ---------------- PANEL: PEAK radius across days (boxplot + swarm) ----------------
peakByDayR = cell(2, numel(BehavCond));
for gi = 1:2
    for d = 1:numel(BehavCond)
        C = radius_bySess{gi, d}; C = C(~cellfun(@isempty, C));
        peakByDayR{gi, d} = cellfun(@(v) max(v(postIdx)), C);
    end
end
allNr=[]; gNr=[]; allTr=[]; gTr=[];
for d = 1:numel(BehavCond)
    allNr=[allNr peakByDayR{1,d}]; gNr=[gNr d*ones(1,numel(peakByDayR{1,d}))]; %#ok<AGROW>
    allTr=[allTr peakByDayR{2,d}]; gTr=[gTr d*ones(1,numel(peakByDayR{2,d}))]; %#ok<AGROW>
end
pKW_Nr = kruskalwallis(allNr, gNr, 'off');
pKW_Tr = kruskalwallis(allTr, gTr, 'off');
fprintf('\n[PEAK radius day effect] Naive KW p=%.4g | Trained KW p=%.4g\n', pKW_Nr, pKW_Tr);
peakByDayPanel(peakByDayR, gngAngles, colNaive, colTrain, ...
    'Peak radius (RMS)', ...
    sprintf('Peak radius by day (KW day effect: Naive p=%.3g, Trained p=%.3g)', pKW_Nr, pKW_Tr), ...
    figDir, 'main_peakRadius_byDay');
try
    fidr = fopen(fullfile(geomTimingDir,'peakRadius_byDay_stats.txt'),'w');
    fprintf(fidr,'Peak radius (max over t>=0, RMS-norm.), per session, by day.\n');
    fprintf(fidr,'Kruskal-Wallis effect of day within group: Naive p=%.4g ; Trained p=%.4g\n\n', pKW_Nr, pKW_Tr);
    fprintf(fidr,'%-8s %-8s %-8s %-8s %-6s\n','group','day','deg','median','n');
    gg={'Naive','Trained'};
    for gi=1:2
        for d=1:numel(BehavCond)
            v=peakByDayR{gi,d};
            fprintf(fidr,'%-8s %-8s %-8d %-8.3f %-6d\n', gg{gi}, BehavCond{d}, gngAngles(d), median(v), numel(v));
        end
    end
    fclose(fidr);
catch
end

%% ---------------- SUSTAINED-WINDOW pooled group comparison (append to stats) ----------------
swMetrics = {'centDist', centDist_bySess, 'peakCentroidDistance_byDay_stats.txt'; ...
             'radius',   radius_bySess,   'peakRadius_byDay_stats.txt'};
swWins = {'1-3s', tvec>=1 & tvec<=3; '0-3s', tvec>=0 & tvec<=3};
for mi = 1:size(swMetrics,1)
    C = swMetrics{mi,2};
    fidS = fopen(fullfile(geomTimingDir, swMetrics{mi,3}), 'a');
    fprintf(fidS, '\n--- Pooled group comparison (all sessions across days), naive vs trained ---\n');
    fprintf(fidS, '    (per-session mean over window; rank-sum two-sided; Hedges g>0 => trained>naive)\n');
    for wi = 1:size(swWins,1)
        widx = swWins{wi,2}; perN=[]; perT=[];
        for d = 1:numel(BehavCond)
            for k=1:numel(C{1,d}), v=C{1,d}{k}; if ~isempty(v), perN(end+1)=mean(v(widx)); end, end %#ok<AGROW>
            for k=1:numel(C{2,d}), v=C{2,d}{k}; if ~isempty(v), perT(end+1)=mean(v(widx)); end, end %#ok<AGROW>
        end
        fprintf(fidS, '  %-5s: Naive %.3f +/- %.3f (n=%d) | Trained %.3f +/- %.3f (n=%d) | rank-sum p=%.4g | Hedges g=%+.2f\n', ...
            swWins{wi,1}, mean(perN), std(perN)/sqrt(numel(perN)), numel(perN), ...
            mean(perT), std(perT)/sqrt(numel(perT)), numel(perT), ranksum(perN,perT), hedgesg(perN,perT));
    end
    fclose(fidS);
end

%% ---------------- LOCAL FUNCTIONS ----------------
function g = hedgesg(x, y)
    x=x(~isnan(x)); y=y(~isnan(y)); n1=numel(x); n2=numel(y);
    sp=sqrt(((n1-1)*var(x)+(n2-1)*var(y))/(n1+n2-2));
    g=((mean(y)-mean(x))/sp)*(1-3/(4*(n1+n2)-9));
end


function peakByDayPanel(peakByDay, gngAngles, colNaive, colTrain, ylab, ttl, figDir, fileBase)
    nD = numel(gngAngles); dayGap = 0.72; off = 0.11; w = 0.14;
    xc = (1:nD)*dayGap;
    fig = figure('Color','w','Position',[100 100 640 440]); ax = axes(fig); hold(ax,'on'); box(ax,'off');
    medN = cellfun(@median, peakByDay(1,:)); medT = cellfun(@median, peakByDay(2,:));
    plot(ax, xc-off, medN, '-', 'Color',[colNaive 0.45], 'LineWidth',1.1);
    plot(ax, xc+off, medT, '-', 'Color',[colTrain 0.45], 'LineWidth',1.1);
    for d = 1:nD
        plotBoxSwarm(ax, xc(d)-off, peakByDay{1,d}, colNaive, w);
        plotBoxSwarm(ax, xc(d)+off, peakByDay{2,d}, colTrain, w);
    end
    set(ax,'XTick',xc,'XTickLabel',arrayfun(@(a)sprintf('%d\\circ',a),gngAngles,'uni',0));
    xlim(ax,[xc(1)-0.40 xc(end)+0.40]);
    xlabel(ax,'Go/NoGo angle'); ylabel(ax,ylab); title(ax, ttl);
    hN = patch(ax,NaN,NaN,colNaive,'FaceAlpha',0.30,'EdgeColor','none');
    hT = patch(ax,NaN,NaN,colTrain,'FaceAlpha',0.30,'EdgeColor','none');
    legend(ax,[hN hT],{'Naive','Trained'},'Box','off','Location','northeast');
    set(ax,'FontSize',11,'LineWidth',1,'TickDir','out');
    save_figure_both(fig, figDir, fileBase);
end

function plotBoxSwarm(ax, x, vals, col, w)
    vals = vals(:); vals = vals(~isnan(vals)); if isempty(vals), return; end
    q = prctile(vals,[25 50 75]); iqr = q(3)-q(1);
    patch(ax, x+[-1 1 1 -1]*w/2, [q(1) q(1) q(3) q(3)], col, 'FaceAlpha',0.30, 'EdgeColor',col, 'LineWidth',0.75);
    line(ax, x+[-1 1]*w/2, [q(2) q(2)], 'Color','k', 'LineWidth',1.1);
    lo = min(vals(vals >= q(1)-1.5*iqr)); hi = max(vals(vals <= q(3)+1.5*iqr));
    line(ax,[x x],[lo q(1)],'Color','k','LineWidth',0.6); line(ax,[x x],[q(3) hi],'Color','k','LineWidth',0.6);
    jx = x + (rand(numel(vals),1)-0.5)*w*0.55;
    scatter(ax, jx, vals, 11, col, 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeColor','none');
end


function save_figure_both(figH, outDir, baseName)
    if nargin < 1 || isempty(figH), figH = gcf; end
    if nargin < 2 || isempty(outDir), outDir = pwd; end
    if nargin < 3 || isempty(baseName), baseName = 'figure'; end
    baseName = regexprep(baseName, '[^\w\-]+', '_');
    figPath = fullfile(outDir, [baseName '.fig']);
    pdfPath = fullfile(outDir, [baseName '.pdf']);
    savefig(figH, figPath);
    try
        exportgraphics(figH, pdfPath, 'ContentType','vector', 'BackgroundColor','white');
    catch
        set(figH, 'PaperPositionMode','auto');
        print(figH, pdfPath, '-dpdf', '-vector');
    end
end

function c = blend_with_white(baseCol, fracWhite)
    fracWhite = min(max(fracWhite,0),1);
    c = (1-fracWhite)*baseCol + fracWhite*[1 1 1];
end

function [dayMean, daySEM] = day_stats(sessCellByDayRow, nDays)
    dayMean = NaN(nDays, 0); daySEM = NaN(nDays, 0);
    nWin = [];
    for d = 1:nDays
        C = sessCellByDayRow{d};
        if isempty(C), continue; end
        C = C(~cellfun(@isempty, C));
        if isempty(C), continue; end
        nWin = numel(C{1}); break;
    end
    if isempty(nWin), return; end
    dayMean = NaN(nDays, nWin); daySEM = NaN(nDays, nWin);
    for d = 1:nDays
        C = sessCellByDayRow{d};
        if isempty(C), continue; end
        C = C(~cellfun(@isempty, C));
        if isempty(C), continue; end
        L = cellfun(@numel, C); C = C(L == nWin);
        if isempty(C), continue; end
        M = vertcat(C{:});
        dayMean(d,:) = mean(M, 1, 'omitnan');
        nEff = sum(isfinite(M), 1);
        sd   = std(M, 0, 1, 'omitnan');
        daySEM(d,:) = sd ./ sqrt(max(nEff,1));
    end
end

function [pooledMean, pooledSEM] = pooled_stats(sessCellByDayRow)
    allSess = {};
    for d = 1:numel(sessCellByDayRow)
        if isempty(sessCellByDayRow{d}), continue; end
        allSess = [allSess, sessCellByDayRow{d}]; %#ok<AGROW>
    end
    allSess = allSess(~cellfun(@isempty, allSess));
    if isempty(allSess), pooledMean = []; pooledSEM = []; return; end
    nWin = numel(allSess{1});
    L = cellfun(@numel, allSess); allSess = allSess(L == nWin);
    if isempty(allSess), pooledMean = []; pooledSEM = []; return; end
    M = vertcat(allSess{:});
    pooledMean = mean(M, 1, 'omitnan');
    nEff = sum(isfinite(M), 1);
    sd   = std(M, 0, 1, 'omitnan');
    pooledSEM = sd ./ sqrt(max(nEff,1));
end

function plot_days_with_alpha(ax, tvec, dayMean, daySEM, lineCol, bandCol, alphaBand, dayAlphas, lw)
% SEM bands drawn first (all days), then mean lines on top.
    nDays = size(dayMean,1);
    % Pass 1: SEM bands
    for d = 1:nDays
        m = dayMean(d,:); sem = daySEM(d,:);
        if all(~isfinite(m)), continue; end
        lo = m - sem; hi = m + sem;
        fill(ax, [tvec fliplr(tvec)], [lo fliplr(hi)], bandCol, ...
            'FaceAlpha', alphaBand, 'EdgeColor','none');
    end
    % Pass 2: mean lines
    for d = 1:nDays
        m = dayMean(d,:);
        if all(~isfinite(m)), continue; end
        plot(ax, tvec, m, '-', 'Color', [lineCol, dayAlphas(d)], 'LineWidth', lw);
    end
end

function plot_mean_with_sem(ax, tvec, m, s, lineCol, bandCol, alphaBand, lw)
    if isempty(m) || isempty(s), return; end
    lo = m - s; hi = m + s;
    fill(ax, [tvec fliplr(tvec)], [lo fliplr(hi)], bandCol, ...
        'FaceAlpha', alphaBand, 'EdgeColor','none');
    plot(ax, tvec, m, '-', 'Color', lineCol, 'LineWidth', lw);
end