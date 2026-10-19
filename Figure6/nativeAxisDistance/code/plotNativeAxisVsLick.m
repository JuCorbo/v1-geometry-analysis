%% plotNativeAxisVsLick.m
%  Figure 6 NATIVE-SPACE supplement. Identical to the PaCMAP axis-position
%  panels (g/h/i) + per-day supp, but the Go->NoGo axis position is measured in
%  the full native neural space (getWindowed 'global'). Reads the native cache.
clear; clc;

%% ---------------- CONFIG ----------------
thisDir = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
anaDir  = fileparts(thisDir);

% native-space cache (matches computeNativeAxisVsLick.m output name)
scheme        = 'global';
snapshotTime  = 1.5;
binSizeTrials = 20; stepSizeTrials = 5;
inMat  = fullfile(anaDir,'cache', sprintf('NativeAxisVsLick_%s_snap%.2fs_bin%d-%d.mat', ...
                  scheme, snapshotTime, binSizeTrials, stepSizeTrials));
tag    = sprintf('native_snap%.1f', snapshotTime);
xLab   = 'Go\rightarrowNoGo axis position (native)';

mainDir = fullfile(anaDir,'panels','main'); if ~exist(mainDir,'dir'), mkdir(mainDir); end
suppDir = fullfile(anaDir,'panels','supp'); if ~exist(suppDir,'dir'), mkdir(suppDir); end
dataDir = fullfile(anaDir, 'data');         if ~exist(dataDir,'dir'),  mkdir(dataDir);  end
statsDir = fullfile(anaDir,'stats');        if ~exist(statsDir,'dir'), mkdir(statsDir); end

goStim        = 45;
exampleDay    = "D3";                 % 30 deg Go/NoGo angle
exampleSess   = ["F17_42221_", "F10_42321_"];   % two examples (steep, flatter)
exampleCols   = [0.85 0.20 0.20; 0.20 0.40 0.85];
nRankGrid     = 40;
minBins       = 5;                    % min sliding-window points to use a session
rMaxFit       = 0.9;                  % fit slope over rank [0, rMaxFit]
epsStd        = 1e-6;

S = load(inMat); T = S.TsumRaw; TT = S.TallRawTrials;
days = sort(unique(T.Behav_Cond));

%% ---------------- per-session NoGo curves + slopes ----------------
% one row per session: Day, Animal, Experiment, nogoStim, slope, (curve stored in C)
Exp=strings(0); Day=strings(0); Ani=strings(0); NoGoDeg=[]; Slope=[];
C = struct('exp',{},'day',{},'x',{},'y',{},'r',{},'z',{});

for di = 1:numel(days)
    Td = T(T.Behav_Cond==days(di), :);
    nogoStim = pickNoGo(Td, goStim);
    if isnan(nogoStim), continue; end
    Tn = Td(Td.Stim==nogoStim, :);
    exps = sort(unique(Tn.Experiment));
    for e = 1:numel(exps)
        Te = Tn(Tn.Experiment==exps(e), :);
        [x, o] = sort(Te.distMedian); y = Te.lickProb(o);
        ok = isfinite(x)&isfinite(y); x=x(ok); y=y(ok);
        if numel(y) < minBins, continue; end
        n = numel(y); r = (0:n-1)'/(n-1);
        z = (y - mean(y)) ./ max(std(y), epsStd);
        m = r <= rMaxFit;
        p = polyfit(r(m), z(m), 1);

        Exp(end+1,1)=exps(e); Day(end+1,1)=days(di); %#ok<AGROW>
        Ani(end+1,1)=string(Te.Animal(1)); NoGoDeg(end+1,1)=nogoStim; %#ok<AGROW>
        Slope(end+1,1)=p(1); %#ok<AGROW>
        C(end+1)=struct('exp',exps(e),'day',days(di),'x',x,'y',y,'r',r,'z',z); %#ok<AGROW>
    end
end
Tslope = table(Day, Ani, Exp, NoGoDeg, Slope, 'VariableNames', ...
    {'Day','Animal','Experiment','NoGoDeg','Slope'});
writetable(Tslope, fullfile(dataDir, sprintf('perSession_slopes_%s.csv', tag)));

nNeg = sum(Tslope.Slope<0); nTot = height(Tslope);
pBino = binocdf(nNeg-1, nTot, 0.5, 'upper');

%% ---- permutation null: shuffle Response within each session ----
rng(7); Nperm = 2000;
Vc = cell(nTot,1); Rc = cell(nTot,1);
for i = 1:nTot
    TTe = TT(TT.Experiment==Tslope.Experiment(i) & TT.Visual_Stim~=goStim, :);
    v = TTe.Metric; r = double(TTe.Response); ok = isfinite(v)&isfinite(r);
    Vc{i} = v(ok); Rc{i} = r(ok);
end
nullS = nan(nTot, Nperm);
for p = 1:Nperm
    for i = 1:nTot
        nullS(i,p) = slopeFromTrials(Vc{i}, Rc{i}(randperm(numel(Rc{i}))), ...
            binSizeTrials, stepSizeTrials, minBins, rMaxFit);
    end
end
medObs   = median(Tslope.Slope,'omitnan');
medNull  = median(nullS,1,'omitnan');
pPerm    = (1 + sum(medNull <= medObs)) / (Nperm + 1);
nullBand = prctile(medNull, [2.5 97.5]);   % null band for the GROUP MEDIAN (narrow)
pFloor   = 1/(Nperm+1);
if pPerm <= pFloor + eps, pPermStr = sprintf('<%.0e', pFloor); else, pPermStr = sprintf('=%.1g', pPerm); end
fprintf('Per-session slopes: %d/%d neg (binom p=%.2g) | median=%.3f | perm p%s | null95=[%.2f %.2f]\n', ...
    nNeg, nTot, pBino, medObs, pPermStr, nullBand(1), nullBand(2));

% --- reported statistics (Figure 6 supp: native-space distance vs licking) ---
fid = fopen(fullfile(statsDir,sprintf('nativeAxisVsLick_stats_%s.txt',tag)),'w');
fprintf(fid,'Figure 6 supp - NoGo axis position vs licking, NATIVE space (no dimensionality reduction)\n');
fprintf(fid,'Metric: per-session slope of z(P(lick)) vs near-to-far rank along the Go->NoGo centroid axis;\n');
fprintf(fid,'snapshot nearest %.1f s, sliding bins of %d trials stepped by %d.\n\n', snapshotTime, binSizeTrials, stepSizeTrials);
fprintf(fid,'n sessions            : %d\n', nTot);
fprintf(fid,'median slope          : %.3f\n', medObs);
fprintf(fid,'negative slopes       : %d/%d\n', nNeg, nTot);
fprintf(fid,'binomial test         : p = %.3g (one-sided, excess of negative slopes) | two-sided p = %.3g\n', pBino, min(1,2*pBino));
fprintf(fid,'median vs null        : perm p %s (observed median vs %d-permutation null of the median;\n', pPermStr, Nperm);
fprintf(fid,'                        null shuffles lick outcomes within each session; floor = 1/(Nperm+1) = %.1e)\n', pFloor);
fprintf(fid,'null 95%% band (median): [%.3f, %.3f]\n', nullBand(1), nullBand(2));
fclose(fid);
fprintf('wrote stats to %s\n', fullfile(statsDir,sprintf('nativeAxisVsLick_stats_%s.txt',tag)));

%% ================= MAIN FIGURE: d | e | f =================
fig = figure('Color','w','Position',[60 120 1550 440]);
allExp = string({C.exp}); allDay = string({C.day});
exNoGo = pickNoGo(T(T.Behav_Cond==exampleDay,:), goStim);

% ---- d) two example sessions: P(lick) curve (raw axis) + short density strip ----
wCol = 0.135; xCol = [0.055, 0.225];
for k = 1:numel(exampleSess)
    ci = find(allExp==exampleSess(k) & allDay==exampleDay, 1);
    if isempty(ci), continue; end
    TTe = TT(TT.Experiment==exampleSess(k) & TT.Visual_Stim~=goStim, :);
    v = TTe.Metric; v = v(isfinite(v)); sep = mean(v);

    axC = axes('Position',[xCol(k) 0.40 wCol 0.50]); hold(axC,'on'); box(axC,'off');
    plot(axC, C(ci).x, C(ci).y, '-o', 'Color',exampleCols(k,:), ...
        'MarkerFaceColor',exampleCols(k,:), 'MarkerSize',4, 'LineWidth',1.8);
    xline(axC,0,'--','Color',[0.55 0.55 0.55]); xline(axC,sep,'-','Color',[0.55 0.55 0.55]);
    ylim(axC,[0 1]); axC.XTickLabel=[];
    if k==1, ylabel(axC,'P(lick) NoGo'); else, axC.YTickLabel=[]; end
    text(axC, 0.04, 0.96, strrep(char(exampleSess(k)),'_','\_'), 'Units','normalized', ...
        'VerticalAlignment','top','Interpreter','tex','FontSize',8,'Color',exampleCols(k,:));

    axD = axes('Position',[xCol(k) 0.255 wCol 0.105]); hold(axD,'on'); box(axD,'off');
    [f,xi]=ksdensity(v);
    fill(axD,[xi fliplr(xi)],[f zeros(size(f))], exampleCols(k,:),'FaceAlpha',0.22,'EdgeColor',exampleCols(k,:),'EdgeAlpha',0.5);
    yr = max(f)*0.10;
    plot(axD,[v v]',[zeros(numel(v),1) -yr*ones(numel(v),1)]','-','Color',[0 0 0 0.25],'LineWidth',0.5);
    xline(axD,0,'--','Color',[0.55 0.55 0.55]); xline(axD,sep,'-','Color',[0.55 0.55 0.55]);
    axD.YTick=[]; if k==1, ylabel(axD,'density','FontSize',8); end
    xlabel(axD,'Go\rightarrowNoGo axis (raw)','FontSize',8);
    linkaxes([axC axD],'x'); xlim(axC,[min(v)-0.5 max(v)+0.5]);
    if k==1
        text(axC, 0, 1.06, sprintf('d   %s examples (%g\\circ)', exampleDay, abs(exNoGo-goStim)), ...
            'Units','normalized','FontWeight','bold','VerticalAlignment','bottom','Interpreter','tex');
    end
end

% ---- e) D3 aggregate: Z(P(lick)) vs near->far rank ----
axe = axes('Position',[0.45 0.255 0.22 0.645]); hold(axe,'on'); box(axe,'off');
rGrid = linspace(0,1,nRankGrid);
Cd = C(allDay==exampleDay);
Z = nan(numel(Cd), nRankGrid);
for i = 1:numel(Cd)
    zi = interp1(Cd(i).r, Cd(i).z, rGrid, 'linear', nan);
    Z(i,:) = zi;
    plot(axe, rGrid, zi, '-', 'Color', 0.82*[1 1 1], 'LineWidth',1);
end
CdExp = string({Cd.exp});
for k = 1:numel(exampleSess)
    ei = find(CdExp==exampleSess(k),1);
    if ~isempty(ei)
        plot(axe, rGrid, interp1(Cd(ei).r,Cd(ei).z,rGrid,'linear',nan), '-', ...
            'Color',exampleCols(k,:),'LineWidth',2.5);
        % dashed linear fit (the slope that enters panel f), over rank [0, rMaxFit]
        rr = Cd(ei).r; zz = Cd(ei).z; mfit = rr<=rMaxFit;
        pk = polyfit(rr(mfit), zz(mfit), 1);
        plot(axe, [0 rMaxFit], polyval(pk,[0 rMaxFit]), '--', ...
            'Color',exampleCols(k,:), 'LineWidth',1.5);
    end
end
mu = mean(Z,1,'omitnan'); se = std(Z,0,1,'omitnan')./sqrt(max(sum(isfinite(Z),1),1));
ok = isfinite(mu)&isfinite(se);
fill(axe,[rGrid(ok) fliplr(rGrid(ok))],[mu(ok)-se(ok) fliplr(mu(ok)+se(ok))],[0 0 0], ...
    'EdgeColor','none','FaceAlpha',0.15);
plot(axe, rGrid, mu, 'k-', 'LineWidth',2.5);
yline(axe,0,':k');
xlabel(axe,'Near\rightarrowfar rank'); ylabel(axe,'Z(P(lick)) NoGo');
ylim(axe,[-2 2]);
title(axe, sprintf('e  %s aggregate (n=%d)', exampleDay, numel(Cd)), 'Interpreter','none');

% ---- f) per-session slopes by day (data) + separate group-median test ----
axf = axes('Position',[0.72 0.255 0.17 0.645]); hold(axf,'on'); box(axf,'off');
dcol = lines(numel(days)); nD = numel(days); xGrp = nD+1; xl = [0.5 xGrp+0.6];
yline(axf, 0, 'k--','LineWidth',1);
% individual session slopes by day (the data; each has a wide single-session null)
[~,loc] = ismember(Tslope.Day, days);
for di = 1:nD
    idx = loc==di; if ~any(idx), continue; end
    jit = 0.10*(rand(sum(idx),1)-0.5);
    scatter(axf, di+jit, Tslope.Slope(idx), 30, dcol(di,:), 'filled','MarkerFaceAlpha',0.85);
end
% separate "group median" column: null-of-median (violin) vs observed median
[fk,xk] = ksdensity(medNull); fk = 0.32*fk/max(fk);
fill(axf, [xGrp+fk fliplr(xGrp-fk)], [xk fliplr(xk)], 0.75*[1 1 1], 'EdgeColor','none','FaceAlpha',0.8);
plot(axf, xGrp+[-0.32 0.32], median(medNull)*[1 1], '-', 'Color',0.5*[1 1 1]);
plot(axf, xGrp, medObs, 'o', 'MarkerFaceColor',[0.85 0.15 0.15], 'MarkerEdgeColor','k', 'MarkerSize',8);
text(axf, xGrp-0.4, medObs, 'obs.', 'Color',[0.85 0.15 0.15],'FontSize',7, ...
    'HorizontalAlignment','right','VerticalAlignment','middle');
xline(axf, nD+0.5, ':', 'Color',0.7*[1 1 1]);
axf.XTick=[1:nD, xGrp]; axf.XTickLabel=[cellstr(days(:))' {'group'}]; xlim(axf,xl);
xlabel(axf,'Day  |  group'); ylabel(axf,'Slope Z(P(lick)) vs rank');
title(axf, sprintf('f  median=%.2f, perm p%s', medObs, pPermStr), 'Interpreter','none');
subtitle(axf, sprintf('%d/%d neg (binomial p=%.0e); grey = null of median', nNeg, nTot, pBino), 'Interpreter','none');
% marginal histogram of the observed slopes (distribution of the data)
drawnow; pf=axf.Position; hw=min(0.045,pf(3)*0.22); gap=0.006;
axh=axes('Position',[pf(1)+pf(3)+gap pf(2) hw pf(4)]); hold(axh,'on'); box(axh,'off');
histogram(axh, Tslope.Slope, 18, 'Orientation','horizontal','FaceColor',0.25*[1 1 1],'EdgeColor','none');
axh.YLim=axf.YLim; axh.YTick=[]; axh.XTick=[]; yline(axh,0,'k--'); yline(axh,medObs,'-','Color',[0.85 0.15 0.15]);

exportgraphics(fig, fullfile(mainDir,sprintf('nativeAxisVsLick_%s.pdf',tag)),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig, fullfile(mainDir,sprintf('nativeAxisVsLick_%s.png',tag)),'Resolution',200);

%% ================= SUPPLEMENT (2x3): every day =================
figS = figure('Color','w','Position',[100 80 1500 780]);
tlS  = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
sgtitle(tlS,'NoGo licking vs Go\rightarrowNoGo axis position — every day (session curves + mean)','Interpreter','tex');
for di = 1:numel(days)
    ax = nexttile(tlS,di); hold(ax,'on'); box(ax,'off');
    Cd = C(string({C.day})==days(di));
    if isempty(Cd), title(ax,sprintf('%s (no data)',days(di))); continue; end
    nogoStim = pickNoGo(T(T.Behav_Cond==days(di),:), goStim);
    % common x-grid across sessions this day (robust 5-95%)
    xall = sort(vertcat(Cd.x)); lo=xall(max(1,round(0.05*numel(xall)))); hi=xall(min(numel(xall),round(0.95*numel(xall))));
    xg = linspace(lo,hi,60); Yg = nan(numel(Cd),numel(xg));
    for i=1:numel(Cd)
        plot(ax, Cd(i).x, Cd(i).y, '-', 'Color',0.8*[1 1 1],'LineWidth',0.8);
        [xu,ia]=unique(Cd(i).x); Yg(i,:)=interp1(xu,Cd(i).y(ia),xg,'linear',nan);
    end
    mu=mean(Yg,1,'omitnan');
    plot(ax, xg, mu, 'k-','LineWidth',2.5);
    ylim(ax,[0 1]); xlim(ax,[lo hi]);
    xlabel(ax, xLab); ylabel(ax,'P(lick) NoGo');
    title(ax, sprintf('%s | NoGo=%g\\circ (angle %g\\circ) | n=%d', days(di), nogoStim, abs(nogoStim-goStim), numel(Cd)),'Interpreter','tex');
end
exportgraphics(figS, fullfile(suppDir,sprintf('nativeAxisVsLick_perDay_%s.pdf',tag)),'ContentType','vector','BackgroundColor','white');
exportgraphics(figS, fullfile(suppDir,sprintf('nativeAxisVsLick_perDay_%s.png',tag)),'Resolution',200);

fprintf('Saved panels to %s (main) / %s (supp)\n', mainDir, suppDir);

%% ================= helpers =================
function nogoStim = pickNoGo(Td, goStim)
    nogoStim = NaN;
    stims = unique(Td.Stim); stims = stims(isfinite(stims));
    cand = stims(stims~=goStim);
    if isempty(cand), return; end
    c = arrayfun(@(s) sum(Td.Stim==s), cand);
    [~,im] = max(c); nogoStim = cand(im);
end

function s = slopeFromTrials(v, resp, binN, binStep, minBins, rMax)
    % sort by metric -> sliding-window P(lick) -> rank -> z -> slope over [0,rMax]
    s = NaN; [~,ord] = sort(v,'ascend'); n = numel(v);
    if n < binN, return; end
    lp = [];
    for i1 = 1:binStep:(n-binN+1), lp(end+1,1) = mean(resp(ord(i1:i1+binN-1))==1); end %#ok<AGROW>
    if numel(lp) < minBins, return; end
    r = (0:numel(lp)-1)'/(numel(lp)-1);
    z = (lp-mean(lp))/max(std(lp),1e-9);
    m = r <= rMax; if sum(m)<3, return; end
    p = polyfit(r(m), z(m), 1); s = p(1);
end
