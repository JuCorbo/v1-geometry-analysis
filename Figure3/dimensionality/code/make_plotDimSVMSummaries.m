% =========================================================================
% make_plotDimSVMSummaries.m
%
% Summary figure for neural dimensionality - D_PR only.
% Layout: 1 x 3
%   Panel 1 - Eigenvalue spectra for eigSpecDay (Go+NoGo averaged per
%             session). Faint lines = individual sessions, thick = mean.
%             Naive (gold) vs Trained (blue). Log y-axis.
%   Panel 2 - D_PR timecourse for eigSpecDay only. Faint lines =
%             individual sessions from that day, thick = mean +/- SD.
%   Panel 3 - D_PR grand average (all days pooled). Mean +/- SD ribbon
%             only - no individual session lines.
%
% Requires: dimMetrics struct as produced by SVMPCATimeCourse_updated.m
% =========================================================================

%% ---- Config ------------------------------------------------------------
eigSpecDay   = 2;        % day shown in panels 1 and 2
maxPCplot    = 50;       % number of PCs to show in spectra
tailRange    = [10 40];  % PC range of the tail-slope fit (as in the compute)
saveDir      = fullfile(fileparts(fileparts(mfilename('fullpath'))),'panels','main');
figName      = 'main_dimensionality_summary';

% ---- Timing test params (mirrors geometry scripts) ----------------------
nPerms       = 10000;
minSustained = 1;
minHedgesG   = 0.2;
% D_PR is expected to be LOWER in Trained -> tail = 'left'
dprTail      = 'left';

colNaive = [0.85 0.65 0.13];   % gold
colTrain = [0.22 0.37 0.65];   % blue
alphaInd = 0.12;                % alpha for individual session lines

%% ---- Setup -------------------------------------------------------------
if ~exist(saveDir, 'dir'), mkdir(saveDir); end

dm_sz_raw = size(dimMetrics);
dm_sz     = ones(1, 4);
dm_sz(1:numel(dm_sz_raw)) = dm_sz_raw;

nDays = dm_sz(2);
nDB   = dm_sz(1);

%% ---- Build data --------------------------------------------------------

% --- Panel 1: eigenspectra for eigSpecDay --------------------------------
% Average Go and NoGo latent vectors per session x rep, then collect.
% Also collect the STORED dimensionality metrics (mean of Go and NoGo) at the
% same snapshot, plus the equal-variance flat-spectrum level h = sum(lam^2)/sum(lam)
% (the flat spectrum of width D_PR and height h has the same total variance and
% the same participation ratio as the data) for the panel annotations.
eigData = cell(nDB, 1);   % each cell: [nSpectra x maxPCplot]
annData = struct('pr', [], 'h', [], 'k80', [], 'tail', []);
annData(nDB) = annData;

for dbi = 1:nDB
    spectra = [];
    for sess = 1:dm_sz(3)
        for rep = 1:dm_sz(4)
            dm = dimMetrics(dbi, eigSpecDay, sess, rep);
            if isempty(dm) || isempty(dm.go), continue; end

            % Average the latent vectors across Go and NoGo for each window,
            % then take the snapshot nearest to t = 1s (post-stim, representative)
            [~, wSnap] = min(abs(dm.tvec - 1.0));

            g = dm.go(wSnap); ng = dm.nogo(wSnap);
            lamGoF = g.latent(:);  lamGoF = lamGoF(lamGoF > 0);
            lamNoF = ng.latent(:); lamNoF = lamNoF(lamNoF > 0);
            if isempty(lamGoF) || isempty(lamNoF), continue; end

            annData(dbi).pr(end+1,1)   = mean([g.D_PR      ng.D_PR],      'omitnan');
            annData(dbi).k80(end+1,1)  = mean([g.k80       ng.k80],       'omitnan');
            annData(dbi).tail(end+1,1) = mean([g.tailSlope ng.tailSlope], 'omitnan');
            annData(dbi).h(end+1,1)    = mean([sum(lamGoF.^2)/sum(lamGoF), ...
                                               sum(lamNoF.^2)/sum(lamNoF)]);

            lamGo   = g.latent(1:min(end, maxPCplot));
            lamNoGo = ng.latent(1:min(end, maxPCplot));

            % Average go/nogo element-wise; treat NaN as missing
            lamAvg  = mean([lamGo(:)'; lamNoGo(:)'], 1, 'omitnan');
            lamAvg(lamAvg <= 0) = NaN;
            if all(isnan(lamAvg)), continue; end

            spectra(end+1, 1:numel(lamAvg)) = lamAvg; %#ok<AGROW>
        end
    end
    eigData{dbi} = spectra;
end

% --- Panel 2: D_PR timecourse for eigSpecDay only -----------------------
dprDay = cell(nDB, 1);   % each cell: [nSessions x nWin]

for dbi = 1:nDB
    dprMat = [];
    for sess = 1:dm_sz(3)
        for rep = 1:dm_sz(4)
            dm = dimMetrics(dbi, eigSpecDay, sess, rep);
            if isempty(dm) || isempty(dm.go), continue; end
            nW   = numel(dm.go);
            rowDPR = nan(1, nW);
            for w = 1:nW
                rowDPR(w) = mean([dm.go(w).D_PR, dm.nogo(w).D_PR], 'omitnan');
            end
            if all(isnan(rowDPR)), continue; end
            dprMat(end+1, :) = rowDPR; %#ok<AGROW>
        end
    end
    dprDay{dbi} = dprMat;
end

% --- Panel 3 & 4: grand averages across all days (D_PR, k80, tailSlope) -
dprAll   = cell(nDB, 1);
k80All   = cell(nDB, 1);
tailAll  = cell(nDB, 1);

for dbi = 1:nDB
    dprMat  = [];
    k80Mat  = [];
    tailMat = [];
    for day = 1:nDays
        for sess = 1:dm_sz(3)
            for rep = 1:dm_sz(4)
                dm = dimMetrics(dbi, day, sess, rep);
                if isempty(dm) || isempty(dm.go), continue; end
                nW = numel(dm.go);
                rowDPR  = nan(1, nW);
                rowK80  = nan(1, nW);
                rowTail = nan(1, nW);
                for w = 1:nW
                    rowDPR(w)  = mean([dm.go(w).D_PR,       dm.nogo(w).D_PR],       'omitnan');
                    rowK80(w)  = mean([dm.go(w).k80,        dm.nogo(w).k80],        'omitnan');
                    rowTail(w) = mean([dm.go(w).tailSlope,  dm.nogo(w).tailSlope],  'omitnan');
                end
                if all(isnan(rowDPR)), continue; end
                dprMat(end+1,  :) = rowDPR;  %#ok<AGROW>
                k80Mat(end+1,  :) = rowK80;  %#ok<AGROW>
                tailMat(end+1, :) = rowTail; %#ok<AGROW>
            end
        end
    end
    dprAll{dbi}  = dprMat;
    k80All{dbi}  = k80Mat;
    tailAll{dbi} = tailMat;
end

% Use tvec from any valid entry
tvec = [];
for dbi = 1:nDB
    for day = 1:nDays
        for sess = 1:dm_sz(3)
            for rep = 1:dm_sz(4)
                dm = dimMetrics(dbi, day, sess, rep);
                if ~isempty(dm) && ~isempty(dm.go)
                    tvec = dm.tvec;
                    break;
                end
            end
            if ~isempty(tvec), break; end
        end
        if ~isempty(tvec), break; end
    end
    if ~isempty(tvec), break; end
end

if isempty(tvec)
    error('No valid dimMetrics entries found.');
end

%% ---- Group separation test: D_PR -----------------------------------------
% Panel 2: per-day test on eigSpecDay
% Panel 3: pooled test across all days

% Per-day test (used in panel 2)
statsDay = groupSeparationTest( ...
    dprDay{1}, dprDay{2}, tvec, ...
    'tail',         dprTail, ...
    'nPerms',       nPerms, ...
    'minSustained', minSustained, ...
    'minHedgesG',   minHedgesG, ...
    'correction',   'cluster');

fprintf('\n=== D_PR group separation - D%d ===\n', eigSpecDay);
fprintf('  Onset (sustained): %.3f s\n', statsDay.onsetTime);
fprintf('  Sig timepoints:    %d / %d\n', sum(statsDay.sigMask), numel(statsDay.sigMask));
fprintf('  Sig (sustained):   %d / %d\n', sum(statsDay.sigMaskSust), numel(statsDay.sigMaskSust));

% Pooled test across all days (used in panel 3)
statsPooled = groupSeparationTest( ...
    dprAll{1}, dprAll{2}, tvec, ...
    'tail',         dprTail, ...
    'nPerms',       nPerms, ...
    'minSustained', minSustained, ...
    'minHedgesG',   minHedgesG, ...
    'correction',   'cluster');

fprintf('\n=== D_PR group separation - grand average (all days) ===\n');
fprintf('  Onset (sustained): %.3f s\n', statsPooled.onsetTime);
fprintf('  Sig timepoints:    %d / %d\n', sum(statsPooled.sigMask), numel(statsPooled.sigMask));
fprintf('  Sig (sustained):   %d / %d\n', sum(statsPooled.sigMaskSust), numel(statsPooled.sigMaskSust));

%% ---- Figure ------------------------------------------------------------
fig = figure('Color', 'w', 'Position', [80 80 1450 340]);
tl  = tiledlayout(1, 4, 'Padding', 'compact', 'TileSpacing', 'loose');

dbLabels = {'Naive', 'Trained'};
cols     = {colNaive, colTrain};

% ---- Panel 1: Eigenspectra + dimensionality-measure annotations ---------
% D_PR : dashed step = flat spectrum with the same total variance and PR
% k80  : dotted vertical at the 80%-cumulative-variance component
% tail : dash-dot line with the group-mean stored tailSlope (natural-log
%        units), anchored to the mean spectrum at PC 25, over tailRange
nexttile(tl); hold on;

muSpec = cell(nDB, 1);
for dbi = 1:nDB
    col      = cols{dbi};
    spectra  = eigData{dbi};
    if isempty(spectra), continue; end

    for k = 1:size(spectra, 1)
        row = spectra(k, :);
        plot(1:size(spectra,2), row, '-', ...
            'Color', [col alphaInd], 'LineWidth', 0.6, 'HandleVisibility', 'off');
    end
    muSpec{dbi} = mean(spectra, 1, 'omitnan');
    plot(1:numel(muSpec{dbi}), muSpec{dbi}, '-', 'Color', col, 'LineWidth', 2.2, ...
        'DisplayName', dbLabels{dbi});
end

set(gca, 'YScale', 'log');
xlim([1 maxPCplot]);
yl1 = ylim; ylim(yl1);

for dbi = 1:nDB
    if isempty(annData(dbi).pr), continue; end
    colD = cols{dbi} * 0.72;
    pr  = mean(annData(dbi).pr);  h  = mean(annData(dbi).h);
    k80 = mean(annData(dbi).k80); tsl = mean(annData(dbi).tail);
    % D_PR: equal-variance flat spectrum
    plot([1 pr], [h h],       '--', 'Color', colD, 'LineWidth', 1.6, 'HandleVisibility', 'off');
    plot([pr pr], [yl1(1) h], '--', 'Color', colD, 'LineWidth', 1.6, 'HandleVisibility', 'off');
    text(pr - 0.5, h*1.12, sprintf('D_{PR} = %.1f', pr), 'Color', colD, 'FontSize', 8, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    % k80: dotted vertical, staggered labels
    plot([k80 k80], yl1, ':', 'Color', colD, 'LineWidth', 1.4, 'HandleVisibility', 'off');
    kLabY = yl1(2) * [0.30 0.52];
    text(k80 + 0.5, kLabY(dbi), sprintf('k_{80} = %.1f', k80), 'Color', colD, 'FontSize', 8, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    % tail slope: stored mean slope anchored at PC 25
    xf = tailRange(1):tailRange(2);
    yf = exp(log(muSpec{dbi}(25)) + tsl * (xf - 25));
    plot(xf, yf, '-.', 'Color', colD*0.75, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    text(tailRange(2) + 0.7, yf(end), sprintf('tail slope = %.3f', tsl), 'Color', colD, ...
        'FontSize', 7.5, 'VerticalAlignment', 'middle');
end

xlabel('PC #');
ylabel('\lambda');
title(sprintf('Eigenvalue spectrum - D%d (t \\approx 1s)', eigSpecDay));
box off; grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1);

% ---- Panel 2: D_PR for eigSpecDay, with individual session lines --------
nexttile(tl); hold on;
nT = numel(tvec);

for dbi = 1:nDB
    col = cols{dbi};
    mat = dprDay{dbi};
    if isempty(mat), continue; end
    if size(mat,2) > nT, mat = mat(:,1:nT); elseif size(mat,2) < nT, mat(:,end+1:nT) = NaN; end

    for k = 1:size(mat, 1)
        plot(tvec, mat(k,:), '-', 'Color', [col alphaInd], 'LineWidth', 0.5);
    end
    mu = mean(mat, 1, 'omitnan');
    sd = std( mat, 0, 1, 'omitnan');
    fill([tvec fliplr(tvec)], [mu+sd fliplr(mu-sd)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(tvec, mu, '-', 'Color', col, 'LineWidth', 2.2, 'DisplayName', dbLabels{dbi});
end

xline(0,   'k--', 'LineWidth', 0.8);
xline(3.0, 'k--', 'LineWidth', 0.8);
xlim([min(tvec) max(tvec)]);
xlabel('Time from stim onset (s)');
ylabel('D_{PR}');
title(sprintf('D_{PR} - D%d, Go+NoGo avg', eigSpecDay));
box off; set(gca, 'FontSize', 11, 'LineWidth', 1);

% Significance bar at top of panel 2
yl2 = ylim;
sigT2 = statsDay.tvec(statsDay.sigMaskSust);
if ~isempty(sigT2)
    yBar2 = yl2(2) + 0.03 * diff(yl2);
    fill([min(sigT2) max(sigT2) max(sigT2) min(sigT2)], ...
         [yBar2-0.012*diff(yl2) yBar2-0.012*diff(yl2) yBar2+0.012*diff(yl2) yBar2+0.012*diff(yl2)], ...
         'k', 'EdgeColor', 'none', 'HandleVisibility', 'off', 'Clipping', 'off');
end
if isfinite(statsDay.onsetTime)
    xline(statsDay.onsetTime, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.4, ...
        'HandleVisibility', 'off');
end

% ---- Panel 3: D_PR grand average, no individual lines -------------------
nexttile(tl); hold on;

for dbi = 1:nDB
    col = cols{dbi};
    mat = dprAll{dbi};
    if isempty(mat), continue; end
    if size(mat,2) > nT, mat = mat(:,1:nT); elseif size(mat,2) < nT, mat(:,end+1:nT) = NaN; end

    mu = mean(mat, 1, 'omitnan');
    sd = std( mat, 0, 1, 'omitnan');
    fill([tvec fliplr(tvec)], [mu+sd fliplr(mu-sd)], col, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(tvec, mu, '-', 'Color', col, 'LineWidth', 2.2, 'DisplayName', dbLabels{dbi});
end

xline(0,   'k--', 'LineWidth', 0.8);
xline(3.0, 'k--', 'LineWidth', 0.8);
xlim([min(tvec) max(tvec)]);
xlabel('Time from stim onset (s)');
ylabel('D_{PR}');
title('D_{PR} - grand average, all days');
box off; set(gca, 'FontSize', 11, 'LineWidth', 1);

% Significance bar at top of panel 3
yl3 = ylim;
sigT3 = statsPooled.tvec(statsPooled.sigMaskSust);
if ~isempty(sigT3)
    yBar3 = yl3(2) + 0.03 * diff(yl3);
    fill([min(sigT3) max(sigT3) max(sigT3) min(sigT3)], ...
         [yBar3-0.012*diff(yl3) yBar3-0.012*diff(yl3) yBar3+0.012*diff(yl3) yBar3+0.012*diff(yl3)], ...
         'k', 'EdgeColor', 'none', 'HandleVisibility', 'off', 'Clipping', 'off');
end
if isfinite(statsPooled.onsetTime)
    xline(statsPooled.onsetTime, '-', 'Color', [0.8 0.2 0.2], 'LineWidth', 1.4, ...
        'HandleVisibility', 'off');
end

sgtitle('Dimensionality summary', 'FontWeight', 'bold', 'FontSize', 13);

% ---- Panel 4: Trained - Naive difference timecourses (D_PR, k80, tailSlope)
axL = nexttile(tl);
hold(axL, 'on');

metLeft   = {dprAll,  k80All};
labLeft   = {'D_{PR}', 'k_{80}'};
colLeft   = {[0.20 0.50 0.80], [0.85 0.33 0.10]};

nT = numel(tvec);
diffLeft = nan(2, nT);
for mi = 1:2
    matN = metLeft{mi}{1}; matT = metLeft{mi}{2};
    if isempty(matN) || isempty(matT), continue; end
    muN = mean(matN, 1, 'omitnan'); muT = mean(matT, 1, 'omitnan');
    if numel(muN) > nT, muN = muN(1:nT); elseif numel(muN) < nT, muN(end+1:nT) = NaN; end
    if numel(muT) > nT, muT = muT(1:nT); elseif numel(muT) < nT, muT(end+1:nT) = NaN; end
    diffLeft(mi,:) = muT - muN;
    plot(axL, tvec, diffLeft(mi,:), '-', 'Color', colLeft{mi}, 'LineWidth', 2.0);
end

yline(axL, 0, 'k:', 'LineWidth', 1.0);
xline(axL, 0,   'k--', 'LineWidth', 0.8);
xline(axL, 3.0, 'k--', 'LineWidth', 0.8);
xlim(axL, [min(tvec) max(tvec)]);
xlabel(axL, 'Time from stim onset (s)');
ylabel(axL, 'Trained - Naive (D_{PR}, k_{80})');
title(axL, 'Difference (Trained - Naive)');
box(axL, 'off'); set(axL, 'FontSize', 11, 'LineWidth', 1);

% Right y-axis for tail slope
axR = axes('Position', axL.Position, 'YAxisLocation', 'right', ...
    'Color', 'none', 'XTick', [], 'FontSize', 11, 'LineWidth', 1);
hold(axR, 'on'); box(axR, 'off');
xlim(axR, [min(tvec) max(tvec)]);

colTail = [0.47 0.67 0.19];
matN = tailAll{1}; matT = tailAll{2};
if ~isempty(matN) && ~isempty(matT)
    muN = mean(matN, 1, 'omitnan'); muT = mean(matT, 1, 'omitnan');
    if numel(muN) > nT, muN = muN(1:nT); elseif numel(muN) < nT, muN(end+1:nT) = NaN; end
    if numel(muT) > nT, muT = muT(1:nT); elseif numel(muT) < nT, muT(end+1:nT) = NaN; end
    diffTail = muT - muN;
    plot(axR, tvec, diffTail, '-', 'Color', colTail, 'LineWidth', 2.0);
end
ylabel(axR, 'Trained - Naive (tail slope)', 'Color', colTail);
axR.YColor = colTail;

% --- align left and right y-axes so 0 sits at the same height on both ---
dL = diffLeft(:); dL = dL(isfinite(dL));
if exist('diffTail','var'), dR = diffTail(:); dR = dR(isfinite(dR)); else, dR = 0; end
loL = min([dL(:);0]); hiL = max([dL(:);0]);
loR = min([dR(:);0]); hiR = max([dR(:);0]);
padL = 0.08*(hiL-loL+eps); padR = 0.08*(hiR-loR+eps);
loL = loL-padL; hiL = hiL+padL; loR = loR-padR; hiR = hiR+padR;
f = -loL/(hiL-loL); f = min(max(f,0.05),0.95);   % zero fraction from the two-sided (D_PR/k80) axis
rangeL = max(hiL/(1-f), -loL/f); ylim(axL, [-f*rangeL, (1-f)*rangeL]);
rangeR = max(hiR/(1-f), -loR/f); ylim(axR, [-f*rangeR, (1-f)*rangeR]);

% Inline labels at right edge
xl4 = xlim(axL);
for mi = 1:2
    lastVal = diffLeft(mi, find(isfinite(diffLeft(mi,:)), 1, 'last'));
    if ~isempty(lastVal)
        text(axL, xl4(2)*1.01, lastVal, labLeft{mi}, ...
            'Color', colLeft{mi}, 'FontSize', 9, ...
            'VerticalAlignment', 'middle', 'Clipping', 'off');
    end
end
if exist('diffTail','var')
    lastTail = diffTail(find(isfinite(diffTail), 1, 'last'));
    if ~isempty(lastTail)
        text(axR, xl4(2)*1.01, lastTail, 'tail slope', ...
            'Color', colTail, 'FontSize', 9, ...
            'VerticalAlignment', 'middle', 'Clipping', 'off');
    end
end

%% ---- Save --------------------------------------------------------------
exportgraphics(fig, fullfile(saveDir, [figName '.pdf']), ...
    'ContentType', 'vector', 'BackgroundColor', 'white');
savefig(fig,  fullfile(saveDir, [figName '.fig']));
saveas(fig,   fullfile(saveDir, [figName '.png']));
fprintf('Saved: %s\n', figName);