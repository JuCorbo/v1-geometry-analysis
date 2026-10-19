%% make_capacityFullTimecourse.m
%
% Figure 4 (manifold capacity) MAIN panel: capacity (rolling window) timecourse.
% Twin layout of the SVM figure -- 1x5 (tile 5 stacked): D2 example |
% Naive all-days | Trained all-days | Grand average | (top) group difference
% per day + (bottom) slope distribution.
%
% Data source (cache, DB-free): <sess>_<DB>_capacity_timecourse.mat in
% ../cache/capacityTimecourse (regenerable with pull_rollingCapacity.m) -- each
% holds capacity_rolling (backward window), tvec, meta (meta.BehavCond -> day,
% meta.DB -> group).
%
% Outputs (into ../panels/main, ../data, ../stats):
%   main_capacity_fullData_timecourse.{pdf,png,fig}
%   capacity_fullData_plotdata_long.csv
%   capacity_fullData_stats.txt
%
% Run:  cd(this code dir); make_capacityFullTimecourse

clear; clc;

%% ---- Paths -------------------------------------------------------------
codeDir = fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir = pwd; end
addpath(codeDir);
anaDir  = fileparts(codeDir);
addpath(fullfile(fileparts(fileparts(anaDir)),'_sharedWindowedData','code'));  % canonical groupSeparationTest
panelDir = fullfile(anaDir, 'panels', 'main'); if ~exist(panelDir,'dir'), mkdir(panelDir); end
dataDir  = fullfile(anaDir, 'data');           if ~exist(dataDir,'dir'),  mkdir(dataDir);  end
statsDir = fullfile(anaDir, 'stats');          if ~exist(statsDir,'dir'), mkdir(statsDir); end

cacheDir = fullfile(anaDir, 'cache', 'capacityTimecourse');
assert(isfolder(cacheDir), 'capacityTimecourse cache folder not found.');
fprintf('Capacity cache: %s\n', cacheDir);

dbNaive = 'BatuhanADNaive'; dbTrained = 'BatuhanAngDist';
BehavCond = {'D1','D2','D3','D4','D5','D6'};
nDays = 6; nDB = 2;

%% ---- Config -----------------------------------------------------------
exampleDay = 2;   % D2

colNaiveM   = [0.85 0.65 0.13];
colNaiveL   = [0.85 0.65 0.13 0.20];
colTrainedM = [0.20 0.40 0.80];
colTrainedL = [0.20 0.40 0.80 0.20];
colStim     = [0.4 0.4 0.4];
colOnset    = [0.8 0.2 0.2];
stimOn      = 0;
stimOff     = 3;
dayAlphas   = linspace(1.0, 0.3, 6);

nPerms       = 10000;
minSustained = 2;
minHedgesG   = 0.2;
capTail      = 'right';   % Trained > Naive
onsetOffsetSec = 0.0646;  % dt: report onset at window CENTRE (matches Fig 1 / Fig 3)

%% ---- Load caches, bin by day/group -> allM{dbi,day} = [nSess x nWin] --
fN = dir(fullfile(cacheDir, sprintf('*_%s_capacity_timecourse.mat', dbNaive)));
fT = dir(fullfile(cacheDir, sprintf('*_%s_capacity_timecourse.mat', dbTrained)));
allFiles = [fN; fT];
assert(~isempty(allFiles), 'No plain *_capacity_timecourse.mat caches found.');

tCommon = [];
raw = struct('dbi',{},'day',{},'roll',{},'tv',{});
for F = allFiles'
    L = load(fullfile(F.folder, F.name), 'capacity_rolling', 'meta', 'tvec');
    if ~isfield(L,'meta') || ~isfield(L.meta,'BehavCond'), continue; end
    day = find(strcmp(BehavCond, L.meta.BehavCond), 1);
    if isempty(day), continue; end
    if contains(L.meta.DB, 'Naive'), dbi = 1; else, dbi = 2; end
    if isempty(tCommon), tCommon = L.tvec(:)'; end
    raw(end+1) = struct('dbi',dbi,'day',day,'roll',L.capacity_rolling(:)','tv',L.tvec(:)'); %#ok<AGROW>
end
t    = tCommon(:)';
nWin = numel(t);

allM = cell(nDB, nDays);
for r = 1:numel(raw)
    v = raw(r).roll;
    if numel(raw(r).tv) ~= nWin || any(abs(raw(r).tv - t) > 1e-9)
        v = interp1(raw(r).tv, v, t, 'linear', 'extrap');
    end
    allM{raw(r).dbi, raw(r).day} = [allM{raw(r).dbi, raw(r).day}; v(:)'];
end
fprintf('Loaded %d sessions across %d files.\n', numel(raw), numel(allFiles));

% Per-session 0-3 s mean capacity (for the by-angle box+swarm summary panel)
gngAngles = [90 45 30 25 20 15];
win03 = t >= 0 & t <= 3;
capByDay = cell(nDB, nDays);
for dbi = 1:nDB
    for day = 1:nDays
        M = allM{dbi, day};
        if ~isempty(M), capByDay{dbi, day} = mean(M(:, win03), 2, 'omitnan'); end
    end
end

%% ---- Per-day group separation tests -----------------------------------
statsDay = cell(1, nDays);
for day = 1:nDays
    MN = allM{1, day}; MT = allM{2, day};
    if isempty(MN) || isempty(MT), continue; end
    statsDay{day} = groupSeparationTest(MN, MT, t, ...
        'tail', capTail, 'correction', 'cluster', 'nPerms', nPerms, ...
        'minSustained', minSustained, 'minHedgesG', minHedgesG, ...
        'onsetOffsetSec', onsetOffsetSec);
    fprintf('D%d: onset = %.3f s  |  nSigSust = %d\n', day, ...
        statsDay{day}.onsetTime, sum(statsDay{day}.sigMaskSust));
end

%% ---- Grand average ----------------------------------------------------
grandN = NaN(nDays, nWin); grandT = NaN(nDays, nWin);
for day = 1:nDays
    if ~isempty(allM{1, day}), grandN(day,:) = mean(allM{1, day}, 1, 'omitnan'); end
    if ~isempty(allM{2, day}), grandT(day,:) = mean(allM{2, day}, 1, 'omitnan'); end
end
muN  = mean(grandN, 1, 'omitnan');
semN = std(grandN,  0, 1, 'omitnan') ./ sqrt(sum(~isnan(grandN), 1));
muT  = mean(grandT, 1, 'omitnan');
semT = std(grandT,  0, 1, 'omitnan') ./ sqrt(sum(~isnan(grandT), 1));

%% ---- Pooled test ------------------------------------------------------
MN_pooled = vertcat(allM{1, :});
MT_pooled = vertcat(allM{2, :});
statsPooled = groupSeparationTest(MN_pooled, MT_pooled, t, ...
    'tail', capTail, 'correction', 'cluster', 'nPerms', nPerms, ...
    'minSustained', minSustained, 'minHedgesG', minHedgesG, ...
    'onsetOffsetSec', onsetOffsetSec, 'verbose', true);
fprintf('Pooled: onset = %.3f s\n', statsPooled.onsetTime);

%% ---- Figure: 1 x 6 ----------------------------------------------------
fig = figure('Color', 'w', 'Position', [80 80 2100 330]);
tl  = tiledlayout(1, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
yLimCap    = [0 0.15];   % example panel (keeps individual session lines in view)
yLimCapAvg = [0 0.1];    % mean panels (Naive/Trained all-days, grand average)
% difference-panel y-limits: tight to the per-day difference range (+margin)
ddiff = grandT - grandN; ddiff = ddiff(isfinite(ddiff));
dpad  = 0.12 * (max(ddiff) - min(ddiff));
yLimDiff = [min(ddiff)-dpad, max(ddiff)+dpad];

% ----- Panel 1: Capacity by Go/NoGo angle (box+swarm) -----
axBox = nexttile(tl); hold(axBox,'on'); box(axBox,'off');
drawCapByDay(axBox, capByDay, gngAngles, colNaiveM, colTrainedM, 'Manifold capacity 0-3 s', 'By Go/NoGo angle');

% ----- Panel 2: Example day (D2) -----
nexttile(tl); hold on; box off;
MN_ex = allM{1, exampleDay}; MT_ex = allM{2, exampleDay};
if ~isempty(MN_ex)
    for k = 1:size(MN_ex,1), plot(t, MN_ex(k,:), '-', 'Color', colNaiveL, 'LineWidth', 0.7); end
end
if ~isempty(MT_ex)
    for k = 1:size(MT_ex,1), plot(t, MT_ex(k,:), '-', 'Color', colTrainedL, 'LineWidth', 0.7); end
end
if ~isempty(MN_ex)
    mu = mean(MN_ex,1,'omitnan'); sem = std(MN_ex,0,1,'omitnan')./sqrt(size(MN_ex,1));
    fill([t fliplr(t)], [mu+sem fliplr(mu-sem)], colNaiveM, 'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
    plot(t, mu, '-', 'Color', colNaiveM, 'LineWidth', 2.4, 'DisplayName','Naive');
end
if ~isempty(MT_ex)
    mu = mean(MT_ex,1,'omitnan'); sem = std(MT_ex,0,1,'omitnan')./sqrt(size(MT_ex,1));
    fill([t fliplr(t)], [mu+sem fliplr(mu-sem)], colTrainedM, 'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
    plot(t, mu, '-', 'Color', colTrainedM, 'LineWidth', 2.4, 'DisplayName','Trained');
end
if ~isempty(statsDay{exampleDay})
    s = statsDay{exampleDay};
    draw_sig_bar(gca, s.tvec, s.sigMaskSust, yLimCap(2));
    if isfinite(s.onsetTime)
        xline(s.onsetTime, '-', 'Color', colOnset, 'LineWidth', 1.2, 'HandleVisibility','off');
    end
end
xline(stimOn,  '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xline(stimOff, '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xlim([-1 4]); ylim(yLimCap);
xlabel('Time (s)'); ylabel('Manifold capacity');
title(sprintf('D%d (example)', exampleDay));

% ----- Panel 2: Naive all days -----
nexttile(tl); hold on; box off;
for day = 1:nDays
    M = allM{1, day}; if isempty(M), continue; end
    mu = mean(M,1,'omitnan'); sem = std(M,0,1,'omitnan')./sqrt(size(M,1));
    fill([t fliplr(t)], [mu+sem fliplr(mu-sem)], colNaiveM, 'EdgeColor','none','FaceAlpha',0.20,'HandleVisibility','off');
end
for day = 1:nDays
    M = allM{1, day}; if isempty(M), continue; end
    mu = mean(M,1,'omitnan');
    plot(t, mu, '-', 'Color', [colNaiveM dayAlphas(day)], 'LineWidth', 2.0, 'DisplayName', sprintf('D%d',day));
end
xline(stimOn,  '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xline(stimOff, '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xlim([-1 4]); ylim(yLimCapAvg);
xlabel('Time (s)'); ylabel('Manifold capacity'); title('Naive - all days');

% ----- Panel 3: Trained all days -----
nexttile(tl); hold on; box off;
for day = 1:nDays
    M = allM{2, day}; if isempty(M), continue; end
    mu = mean(M,1,'omitnan'); sem = std(M,0,1,'omitnan')./sqrt(size(M,1));
    fill([t fliplr(t)], [mu+sem fliplr(mu-sem)], colTrainedM, 'EdgeColor','none','FaceAlpha',0.20,'HandleVisibility','off');
end
for day = 1:nDays
    M = allM{2, day}; if isempty(M), continue; end
    mu = mean(M,1,'omitnan');
    plot(t, mu, '-', 'Color', [colTrainedM dayAlphas(day)], 'LineWidth', 2.0, 'DisplayName', sprintf('D%d',day));
end
xline(stimOn,  '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xline(stimOff, '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xlim([-1 4]); ylim(yLimCapAvg);
xlabel('Time (s)'); ylabel('Manifold capacity'); title('Trained - all days');

% ----- Panel 4: Grand average + pooled significance -----
ax4 = nexttile(tl); hold on; box off;
fill([t fliplr(t)], [muN+semN fliplr(muN-semN)], colNaiveM, 'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
fill([t fliplr(t)], [muT+semT fliplr(muT-semT)], colTrainedM, 'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
plot(t, muN, '-', 'Color', colNaiveM,   'LineWidth', 2.4, 'DisplayName','Naive');
plot(t, muT, '-', 'Color', colTrainedM, 'LineWidth', 2.4, 'DisplayName','Trained');
draw_sig_bar(ax4, statsPooled.tvec, statsPooled.sigMaskSust, yLimCapAvg(2));
if isfinite(statsPooled.onsetTime)
    xline(statsPooled.onsetTime, '-', 'Color', colOnset, 'LineWidth', 1.2, 'HandleVisibility','off');
end
xline(stimOn,  '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xline(stimOff, '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xlim([-1 4]); ylim(yLimCapAvg);
xlabel('Time (s)'); ylabel('Manifold capacity'); title('Grand average');

% ----- Panel 5: stacked (difference top + slope histogram bottom) -----
axDummy = nexttile(tl); pos = axDummy.Position; delete(axDummy);
hGap = 0.04*pos(4); avail = pos(4)-hGap; hTop = 0.60*avail; hBot = 0.40*avail; yBot = pos(2); yTop = pos(2)+hBot+hGap;

axDiff = axes(fig, 'Position', [pos(1) yTop pos(3) hTop]); hold(axDiff,'on'); box(axDiff,'off');
yline(axDiff, 0, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
for day = 1:nDays
    mN = grandN(day,:); mT = grandT(day,:);
    if all(~isfinite(mN)) || all(~isfinite(mT)), continue; end
    plot(axDiff, t, mT - mN, '-', 'Color', [0 0 0 dayAlphas(day)], 'LineWidth', 2.0);
end
xline(axDiff, stimOn,  '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xline(axDiff, stimOff, '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
xlim(axDiff, [-1 4]); ylim(axDiff, yLimDiff);
ylabel(axDiff, '\Delta Capacity'); title(axDiff, 'Group difference per day');
set(axDiff, 'XTickLabel', {}, 'FontSize', 10, 'LineWidth', 1);

fitWin = [1, 3];
slopeNaive_cap = []; slopeTrained_cap = [];
for day = 1:nDays
    idx = t >= fitWin(1) & t <= fitWin(2);
    MN = allM{1, day};
    if ~isempty(MN)
        for k = 1:size(MN,1)
            tc = MN(k,:); if sum(isfinite(tc(idx))) < 3, continue; end
            p = polyfit(t(idx), tc(idx), 1); slopeNaive_cap(end+1) = p(1); %#ok<AGROW>
        end
    end
    MT = allM{2, day};
    if ~isempty(MT)
        for k = 1:size(MT,1)
            tc = MT(k,:); if sum(isfinite(tc(idx))) < 3, continue; end
            p = polyfit(t(idx), tc(idx), 1); slopeTrained_cap(end+1) = p(1); %#ok<AGROW>
        end
    end
end
axHist = axes(fig, 'Position', [pos(1) yBot pos(3) hBot]); hold(axHist,'on'); box(axHist,'off');
allSlopes = [slopeNaive_cap(:); slopeTrained_cap(:)];
if ~isempty(allSlopes)
    binEdges = linspace(min(allSlopes), max(allSlopes), 20);
    histogram(axHist, slopeNaive_cap,   binEdges, 'FaceColor', colNaiveM,   'FaceAlpha', 0.6, 'EdgeColor','none');
    histogram(axHist, slopeTrained_cap, binEdges, 'FaceColor', colTrainedM, 'FaceAlpha', 0.6, 'EdgeColor','none');
end
xline(axHist, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, 'HandleVisibility','off');
xlabel(axHist, sprintf('Capacity slope %.0f-%.0fs (/s)', fitWin(1), fitWin(2)));
ylabel(axHist, 'Count'); title(axHist, 'Slope distribution');
set(axHist, 'FontSize', 10, 'LineWidth', 1);

%% ---- Save figure ------------------------------------------------------
base = fullfile(panelDir, 'main_capacity_fullData_timecourse');
exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
exportgraphics(fig, [base '.png'], 'Resolution', 300, 'BackgroundColor','white');
savefig(fig, [base '.fig']);
fprintf('Saved panel: %s\n', base);

%% ---- Left summary panel: capacity by Go/NoGo angle (box+swarm) --------
capByDayPanel(capByDay, gngAngles, colNaiveM, colTrainedM, ...
    'Manifold capacity 0-3 s', 'Capacity by Go/NoGo angle', ...
    panelDir, 'main_capacity_capacityByDay_boxswarm');

%% ---- SUPP: per-day timecourses (2x3) ----------------------------------
suppDir = fullfile(anaDir, 'panels', 'supp'); if ~exist(suppDir,'dir'), mkdir(suppDir); end
figPD = figure('Color','w','Position',[120 120 1100 640]);
tlPD  = tiledlayout(figPD, 2, 3, 'TileSpacing','compact','Padding','compact');
for day = 1:nDays
    ax = nexttile(tlPD); hold(ax,'on'); box(ax,'off');
    MN = allM{1,day}; MT = allM{2,day};
    for k = 1:size(MN,1), plot(ax, t, MN(k,:), '-', 'Color', colNaiveL,   'LineWidth', 0.6, 'HandleVisibility','off'); end
    for k = 1:size(MT,1), plot(ax, t, MT(k,:), '-', 'Color', colTrainedL, 'LineWidth', 0.6, 'HandleVisibility','off'); end
    if ~isempty(MN)
        mu = mean(MN,1,'omitnan'); sem = std(MN,0,1,'omitnan')./sqrt(size(MN,1));
        fill(ax,[t fliplr(t)],[mu+sem fliplr(mu-sem)],colNaiveM,'EdgeColor','none','FaceAlpha',0.20,'HandleVisibility','off');
        plot(ax,t,mu,'-','Color',colNaiveM,'LineWidth',2.0,'DisplayName','Naive');
    end
    if ~isempty(MT)
        mu = mean(MT,1,'omitnan'); sem = std(MT,0,1,'omitnan')./sqrt(size(MT,1));
        fill(ax,[t fliplr(t)],[mu+sem fliplr(mu-sem)],colTrainedM,'EdgeColor','none','FaceAlpha',0.20,'HandleVisibility','off');
        plot(ax,t,mu,'-','Color',colTrainedM,'LineWidth',2.0,'DisplayName','Trained');
    end
    if ~isempty(statsDay{day})
        s = statsDay{day};
        draw_sig_bar(ax, s.tvec, s.sigMaskSust, yLimCap(2));
        if isfinite(s.onsetTime), xline(ax, s.onsetTime, '-', 'Color', colOnset, 'LineWidth', 1.2, 'HandleVisibility','off'); end
    end
    xline(ax, stimOn,  '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(ax, stimOff, '--', 'Color', colStim, 'LineWidth', 0.8, 'HandleVisibility','off');
    xlim(ax,[-1 4]); ylim(ax, yLimCap);
    xlabel(ax,'Time (s)'); ylabel(ax,'Manifold capacity'); title(ax, sprintf('D%d', day));
    if day == 1, legend(ax,'Location','northeast','Box','off'); end
end
baseS = fullfile(suppDir, 'supp_capacity_perDay_timecourse');
exportgraphics(figPD, [baseS '.pdf'], 'ContentType','vector','BackgroundColor','white');
exportgraphics(figPD, [baseS '.png'], 'Resolution',300,'BackgroundColor','white');
savefig(figPD, [baseS '.fig']);
fprintf('Saved supp per-day: %s\n', baseS);

%% ---- Plot-data CSV (long) ---------------------------------------------
csvPath = fullfile(dataDir, 'capacity_fullData_plotdata_long.csv');
fid = fopen(csvPath, 'w');
fprintf(fid, 'group,day,session,time_s,capacity_rolling\n');
grpName = {'Naive','Trained'};
for dbi = 1:nDB
    for day = 1:nDays
        M = allM{dbi, day}; if isempty(M), continue; end
        for k = 1:size(M,1)
            for w = 1:nWin
                if isfinite(M(k,w))
                    fprintf(fid, '%s,D%d,%d,%.4f,%.6f\n', grpName{dbi}, day, k, t(w), M(k,w));
                end
            end
        end
    end
end
fclose(fid);
fprintf('Saved plot data: %s\n', csvPath);

%% ---- Stats file -------------------------------------------------------
txtPath = fullfile(statsDir, 'capacity_fullData_stats.txt');
fid = fopen(txtPath, 'w');
fprintf(fid, 'Manifold capacity (rolling) timecourse - group separation (Trained vs Naive)\n');
fprintf(fid, 'Source caches: cache/capacityTimecourse\n');
fprintf(fid, 'Test: groupSeparationTest (cluster-based permutation, nPerms=%d, tail=%s, |Hedges g|>=%.2f, minSustained=%d)\n', ...
    nPerms, capTail, minHedgesG, minSustained);
fprintf(fid, 'Onset reported at window CENTRE (onsetOffsetSec=%.4f; matches Fig 1 / Fig 3).\n', onsetOffsetSec);
fprintf(fid, 'clusterP = fraction of 10,000 permutations with max cluster mass >= observed (Maris & Oostenveld 2007), floored at 1e-4.\n\n');
fprintf(fid, '%-6s  %-12s  %-10s  %-10s  %-8s\n', 'Day', 'Onset(s)', 'clusterP', 'nSigSust', 'nSess N/T');
for day = 1:nDays
    if isempty(statsDay{day}), continue; end
    s = statsDay{day};
    fprintf(fid, 'D%-5d  %-12s  %-10.4g  %-10d  %d/%d\n', day, num2str(s.onsetTime,'%.3f'), ...
        s.clusterPmin, sum(s.sigMaskSust), size(allM{1,day},1), size(allM{2,day},1));
end
fprintf(fid, '\nPooled (all sessions x days): onset = %.3f s, clusterP = %.4g, nSigSust = %d, nSess N/T = %d/%d\n', ...
    statsPooled.onsetTime, statsPooled.clusterPmin, sum(statsPooled.sigMaskSust), size(MN_pooled,1), size(MT_pooled,1));
fprintf(fid, '\nSlope 1-3s (capacity per s):\n');
fprintf(fid, '  Naive:   n=%d, mean=%.4f, median=%.4f\n', numel(slopeNaive_cap), mean(slopeNaive_cap), median(slopeNaive_cap));
fprintf(fid, '  Trained: n=%d, mean=%.4f, median=%.4f\n', numel(slopeTrained_cap), mean(slopeTrained_cap), median(slopeTrained_cap));
if ~isempty(slopeNaive_cap) && ~isempty(slopeTrained_cap)
    pRS = ranksum(slopeNaive_cap, slopeTrained_cap);
    fprintf(fid, '  Naive vs Trained slope, rank-sum p = %.4g\n', pRS);
end
fclose(fid);
fprintf('Saved stats: %s\n', txtPath);
fprintf('DONE (capacity).\n');

%% ---- Local helpers ----------------------------------------------------
function drawCapByDay(ax, byDay, gngAngles, colNaive, colTrain, ylab, ttl)
    nD = numel(gngAngles); dayGap = 0.72; off = 0.11; w = 0.14;
    xc = (1:nD)*dayGap;
    medN = cellfun(@median, byDay(1,:)); medT = cellfun(@median, byDay(2,:));
    plot(ax, xc-off, medN, '-', 'Color',[colNaive 0.45], 'LineWidth',1.1, 'HandleVisibility','off');
    plot(ax, xc+off, medT, '-', 'Color',[colTrain 0.45], 'LineWidth',1.1, 'HandleVisibility','off');
    for d = 1:nD
        plotBoxSwarm(ax, xc(d)-off, byDay{1,d}, colNaive, w);
        plotBoxSwarm(ax, xc(d)+off, byDay{2,d}, colTrain, w);
    end
    set(ax,'XTick',xc,'XTickLabel',arrayfun(@(a)sprintf('%d\\circ',a),gngAngles,'uni',0));
    xlim(ax,[xc(1)-0.40 xc(end)+0.40]); ylim(ax,[0 0.09]);
    xlabel(ax,'Go/NoGo angle'); ylabel(ax,ylab); title(ax, ttl);
    hN = patch(ax,NaN,NaN,colNaive,'FaceAlpha',0.30,'EdgeColor','none');
    hT = patch(ax,NaN,NaN,colTrain,'FaceAlpha',0.30,'EdgeColor','none');
    legend(ax,[hN hT],{'Naive','Trained'},'Box','off','Location','northeast');
    set(ax,'FontSize',10,'LineWidth',1,'TickDir','out');
end

function capByDayPanel(byDay, gngAngles, colNaive, colTrain, ylab, ttl, figDir, fileBase)
    fig = figure('Color','w','Position',[100 100 560 440]); ax = axes(fig); hold(ax,'on'); box(ax,'off');
    drawCapByDay(ax, byDay, gngAngles, colNaive, colTrain, ylab, ttl);
    set(ax,'FontSize',11);
    base = fullfile(figDir, fileBase);
    exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
    exportgraphics(fig, [base '.png'], 'Resolution',300, 'BackgroundColor','white');
    savefig(fig, [base '.fig']);
    fprintf('Saved panel: %s\n', base);
end

function plotBoxSwarm(ax, x, vals, col, w)
    vals = vals(:); vals = vals(~isnan(vals)); if isempty(vals), return; end
    q = prctile(vals,[25 50 75]); iqr = q(3)-q(1);
    patch(ax, x+[-1 1 1 -1]*w/2, [q(1) q(1) q(3) q(3)], col, 'FaceAlpha',0.30, 'EdgeColor',col, 'LineWidth',0.75);
    line(ax, x+[-1 1]*w/2, [q(2) q(2)], 'Color','k', 'LineWidth',1.1);
    lo = min(vals(vals >= q(1)-1.5*iqr)); hi = max(vals(vals <= q(3)+1.5*iqr));
    line(ax,[x x],[lo q(1)],'Color','k','LineWidth',0.6); line(ax,[x x],[q(3) hi],'Color','k','LineWidth',0.6);
    jx = x + (rand(numel(vals),1)-0.5)*w*0.55;
    scatter(ax, jx, vals, 14, col, 'filled', 'MarkerFaceAlpha',0.7, 'MarkerEdgeColor','none');
end

function draw_sig_bar(ax, tvec_s, sigMask, yTop)
sigT = tvec_s(sigMask);
if isempty(sigT), return; end
yl = ylim(ax); span = diff(yl);
yBar = yTop + 0.03*span; barH = 0.012*span;
fill(ax, [min(sigT) max(sigT) max(sigT) min(sigT)], ...
     [yBar-barH yBar-barH yBar+barH yBar+barH], ...
     'k', 'EdgeColor','none', 'HandleVisibility','off', 'Clipping','off');
end
