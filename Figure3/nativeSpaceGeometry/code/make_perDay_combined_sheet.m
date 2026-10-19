%% make_perDay_combined_sheet.m
%  Combine the six per-day native-geometry supplements into ONE single-sheet PDF:
%  rows = measure (radius, centroid distance, radius Go, radius NoGo,
%  centroid magnitude Go, centroid magnitude NoGo), columns = day D1..D6.
%  Each cell = naive/trained per-session traces (thin) + group means (thick),
%  x-axis from baseline (~ -1 s), one shared y-scale per row (per measure).
%  Reads the same cache as make_originalSpaceGeometryPlotting.m (SQL-free);
%  helper copies below are byte-identical to that script.
clear; clc;
codeDir = fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir = pwd; end
load(fullfile(fileparts(codeDir),'cache','originalSpaceBasicGeometrySession.mat'));  % data cells, tvec, colours, styles, BehavCond
suppDir = fullfile(fileparts(codeDir),'panels','supp'); if ~exist(suppDir,'dir'), mkdir(suppDir); end

metrics = { ...
    'Radius',                    radius_bySess; ...
    'Centroid distance',         centDist_bySess; ...
    'Radius (Go)',               radiusGo_bySess; ...
    'Radius (NoGo)',             radiusNo_bySess; ...
    'Centroid mag. (Go)',        centMagGo_bySess; ...
    'Centroid mag. (NoGo)',      centMagNo_bySess };
nM = size(metrics,1); nD = numel(BehavCond);

fig = figure('Color','w','Position',[40 20 1700 1650]);
tl  = tiledlayout(nM, nD, 'TileSpacing','compact','Padding','compact');
sgtitle(tl, 'Native-space geometry, per day (rows = measure, columns = day)', 'FontWeight','bold');

for r = 1:nM
    C = metrics{r,2};
    % per-row (per-measure) shared y-limits over all days and both groups
    allY = [];
    for g = 1:2, for d = 1:nD
        cc = C{g,d};
        for i = 1:numel(cc), allY = [allY; cc{i}(:)]; end %#ok<AGROW>
    end, end
    yl = [min(allY) max(allY)];
    if isempty(allY) || ~all(isfinite(yl)) || yl(1)==yl(2), yl = [0 1]; end
    pad = 0.05*(yl(2)-yl(1)); yl = [yl(1)-pad yl(2)+pad];

    for d = 1:nD
        ax = nexttile(tl, (r-1)*nD + d); hold(ax,'on'); box(ax,'off');
        plot_group_sessions(ax, C{1,d}, tvec, colNaive, thinAlpha, thinLW, '-');
        plot_group_sessions(ax, C{2,d}, tvec, colTrain, thinAlpha, thinLW, '-');
        [mN, okN] = group_mean(C{1,d}');
        [mT, okT] = group_mean(C{2,d}');
        if okN, plot(ax, tvec, mN, '-', 'Color', colNaive, 'LineWidth', thickLW); end
        if okT, plot(ax, tvec, mT, '-', 'Color', colTrain, 'LineWidth', thickLW); end
        xline(ax, 0, ':', 'Color', [0.6 0.6 0.6]);
        xlim(ax, [tvec(1) tvec(end)]); ylim(ax, yl);
        set(ax, 'FontSize', 7);
        if r == 1, title(ax, BehavCond{d}, 'Interpreter','none'); end
        if d == 1, ylabel(ax, metrics{r,1}); else, set(ax,'YTickLabel',{}); end
        if r <  nM, set(ax,'XTickLabel',{}); else, xlabel(ax,'Time (s)'); end
        % legend once, in the top-right cell
        if r == 1 && d == nD
            cNt = blend_with_white(colNaive, 1-thinAlpha);
            cTt = blend_with_white(colTrain, 1-thinAlpha);
            h1 = plot(ax, nan, nan, '-', 'Color', cNt,      'LineWidth', thinLW);
            h2 = plot(ax, nan, nan, '-', 'Color', cTt,      'LineWidth', thinLW);
            h3 = plot(ax, nan, nan, '-', 'Color', colNaive, 'LineWidth', thickLW);
            h4 = plot(ax, nan, nan, '-', 'Color', colTrain, 'LineWidth', thickLW);
            legend(ax, [h1 h2 h3 h4], {'naive sessions','trained sessions','naive mean','trained mean'}, ...
                'FontSize', 6, 'Box','off', 'Location','northeast');
        end
    end
end

base = fullfile(suppDir, 'supp_perDay_ALL');
exportgraphics(fig, [base '.pdf'], 'ContentType','vector', 'BackgroundColor','white');
exportgraphics(fig, [base '.png'], 'Resolution', 160, 'BackgroundColor','white');
savefig(fig, [base '.fig']);
fprintf('Saved combined per-day sheet to %s.(pdf/png/fig)\n', base);

%% ---- helpers (copied verbatim from make_originalSpaceGeometryPlotting.m) ----
function plot_group_sessions(ax, sessCell, tvec, baseCol, alphaLike, lw, ls)
    if isempty(sessCell), return; end
    if nargin < 7 || isempty(ls), ls = '-'; end
    colThin = blend_with_white(baseCol, 1-alphaLike);
    for i = 1:numel(sessCell)
        y = sessCell{i};
        if isempty(y) || numel(y) ~= numel(tvec), continue; end
        plot(ax, tvec, y, ls, 'Color', colThin, 'LineWidth', lw);
    end
end

function [m, ok] = group_mean(sessCell)
    ok = false; m = [];
    if isempty(sessCell), return; end
    L = cellfun(@numel, sessCell);
    if any(L == 0) || numel(unique(L)) ~= 1, return; end
    M = cell2mat(cellfun(@(v) v(:)', sessCell, 'UniformOutput', false));
    m = mean(M, 1, 'omitnan');
    ok = true;
end

function c = blend_with_white(baseCol, fracWhite)
    fracWhite = min(max(fracWhite,0),1);
    c = (1-fracWhite)*baseCol + fracWhite*[1 1 1];
end
