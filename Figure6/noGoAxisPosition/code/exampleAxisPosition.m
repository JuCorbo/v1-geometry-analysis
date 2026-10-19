%% exampleAxisPosition.m
%  Figure 6 — panel f. Example single-session PaCMAP clouds; NoGo trials colored by
%  their position on the Go->NoGo centroid axis (Near -> Far, rank-mapped). The
%  Go/NoGo separating hyperplane is drawn for illustration.
%
%  ORIENTATION RULE: PaCMAP axis orientation is arbitrary
%  per run, so each session's snapshot is Procrustes-aligned (rigid: no scaling,
%  best reflection) to a single global REFERENCE session's Go/NoGo centroid-
%  trajectory landmarks, then rendered with the default view(3). This is the
%  codebase's "aligned" rule; it is display-only and rotation-invariant so it does
%  not change any distance.
%
%  Reads the main global-z (seed0) embeddings at the 1.5 s snapshot.
clear; clc;

%% ---------------- CONFIG ----------------
seed           = 0;
snapshotTime   = 1.5;
goVal          = 45;
dimsToUse      = [1 2 3];
sessionsToPlot = ["F17_42221_", "F10_42321_"];     % same examples as panels d/e (D3, 30 deg)
refSession     = "";                                % "" = first valid session as global ref
keepFrac       = 0.90;
cmapName       = 'parula';

thisDir  = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
anaDir   = fileparts(thisDir);
finalRoot= fileparts(fileparts(anaDir));
embRoot  = fullfile(finalRoot,'Figure2','pacmap','cache','main_pacmap_bwgz_pm080', ...
                    sprintf('baseline_seed%d',seed),'BatuhanAngDist');
panelDir = fullfile(anaDir,'panels','main'); if ~exist(panelDir,'dir'), mkdir(panelDir); end

files = [];
for d = 1:6
    files = [files; dir(fullfile(embRoot, sprintf('D%d',d), '*_PaCMAP_UMATO.mat'))]; %#ok<AGROW>
end
names = string({files.name});

%% ---------------- GLOBAL REFERENCE (for orientation) ----------------
ref = buildReference(files, names, refSession, dimsToUse, goVal);
fprintf('Orientation reference: %s\n', ref.name);

%% ---------------- PLOT ----------------
fig = figure('Color','w','Position',[100 100 1400 640]);
tl  = tiledlayout(1, numel(sessionsToPlot), 'Padding','compact','TileSpacing','compact');

for si = 1:numel(sessionsToPlot)
    hit = find(contains(names, sessionsToPlot(si)), 1, 'first');
    assert(~isempty(hit), 'Session %s not found under %s', sessionsToPlot(si), embRoot);
    S = load(fullfile(files(hit).folder, files(hit).name));

    Y   = double(S.Y_pacmap(:, dimsToUse));
    lab = double(S.labels(:));
    twi = double(S.time_window_idx(:));
    twc = double(S.time_window_center(:));
    stimVals = unique(lab); stimVals = stimVals(isfinite(stimVals));
    nogoVal  = stimVals(stimVals ~= goVal); nogoVal = nogoVal(1);

    % --- Procrustes-align this session to the reference frame ---
    Cs = makeLandmarks(Y, lab, twi, twc, goVal, nogoVal);
    good = all(isfinite(ref.Cstack),2) & all(isfinite(Cs),2);
    [~,~,tr] = procrustes(ref.Cstack(good,:), Cs(good,:), 'Scaling',false, 'Reflection','best');
    Ya = tr.b * Y * tr.T + tr.c(1,:);       % aligned coordinates (rigid)

    % --- snapshot window ---
    uW = sort(unique(twi(isfinite(twi))));
    tC = arrayfun(@(w) mean(twc(twi==w),'omitnan'), uW);
    [~,k] = min(abs(tC - snapshotTime)); snapWIdx = uW(k); snapT = tC(k);
    idx = (twi==snapWIdx);
    X = Ya(idx,:); L = lab(idx);

    % --- SVM Go vs NoGo in aligned space (plane = illustration only) ---
    y = nan(size(L)); y(L==goVal)=+1; y(L==nogoVal)=-1;
    mdl = fitcsvm(X, y, 'KernelFunction','linear', 'Standardize',false, 'ClassNames',[-1 +1]);
    w = mdl.Beta; b = mdl.Bias;

    % --- NoGo colored by Go->NoGo axis position, RANK-mapped (density-equalized) ---
    % Raw axisProj is compressed within the cloud (most trials near the NoGo
    % centroid), so a linear colormap would make the bulk one color. Coloring by
    % the quantile of axisProj spreads Near->Far evenly across the visible cloud,
    % consistent with the rank-based analysis.
    idxGo = (L==goVal); idxNo = (L==nogoVal);
    goC = mean(X(idxGo,:),1); nogoC = mean(X(idxNo,:),1);
    u = (nogoC - goC); u = u/norm(u);
    Xno = X(idxNo,:); axNo = (Xno - goC) * u';       % Near (Go-like) -> Far
    qNo = tiedrank(axNo) / numel(axNo);              % density-equalized [0,1]
    if keepFrac < 1                                   % trim the far tail for a tight view
        [~,ord] = sort(axNo,'ascend'); nKeep = max(1, floor(keepFrac*numel(axNo)));
        m = false(size(axNo)); m(ord(1:nKeep)) = true; XnoP = Xno(m,:); corrP = qNo(m);
    else
        XnoP = Xno; corrP = qNo;
    end

    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); view(ax,3);
    xlabel(ax,'PaCMAP1'); ylabel(ax,'PaCMAP2'); zlabel(ax,'PaCMAP3');

    scatter3(ax, X(idxGo,1), X(idxGo,2), X(idxGo,3), 16, 0.35*[1 1 1], 'filled', 'MarkerFaceAlpha',0.5);
    scatter3(ax, XnoP(:,1), XnoP(:,2), XnoP(:,3), 30, corrP, 'filled', 'Marker','^','MarkerFaceAlpha',0.9);

    % SVM plane
    mins = min(X,[],1); maxs = max(X,[],1); pad = 0.1*(maxs-mins); mins=mins-pad; maxs=maxs+pad;
    [~,solveDim] = max(abs(w)); dd = 1:3; dd(dd==solveDim)=[]; d1=dd(1); d2=dd(2); d3=solveDim;
    [U,V] = meshgrid(linspace(mins(d1),maxs(d1),25), linspace(mins(d2),maxs(d2),25));
    Z = -(w(d1)*U + w(d2)*V + b)/w(d3);
    P = zeros(size(U,1),size(U,2),3); P(:,:,d1)=U; P(:,:,d2)=V; P(:,:,d3)=Z;
    surf(ax, P(:,:,1), P(:,:,2), P(:,:,3), 'FaceAlpha',0.12,'EdgeAlpha',0.1,'FaceColor',[0 0 0]);

    xlim(ax, prctile(X(:,1),[5 95]) + [-.05 .05]*range(prctile(X(:,1),[5 95])));
    ylim(ax, prctile(X(:,2),[5 95]) + [-.05 .05]*range(prctile(X(:,2),[5 95])));
    zlim(ax, prctile(X(:,3),[5 95]) + [-.05 .05]*range(prctile(X(:,3),[5 95])));
    axis(ax,'vis3d');

    colormap(ax, feval(cmapName,256));
    cl = [min(corrP) max(corrP)]; if ~all(isfinite(cl))||diff(cl)<=0, cl=[0 1]; end
    caxis(ax,cl); cb = colorbar(ax); cb.Ticks=cl; cb.TickLabels={'Near','Far'};

    title(ax, sprintf('%s   snapshot %.2f s', strrep(char(sessionsToPlot(si)),'_','\_'), snapT), 'Interpreter','tex');
end

exportgraphics(fig, fullfile(panelDir,'fig6_f_exampleAxisPosition.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig, fullfile(panelDir,'fig6_f_exampleAxisPosition.png'),'Resolution',200);
fprintf('Saved panel c to %s\n', panelDir);

%% ================= helpers =================
function ref = buildReference(files, names, refSession, dimsToUse, goVal)
    order = 1:numel(files);
    if strlength(refSession) > 0
        h = find(contains(names, refSession),1,'first');
        if ~isempty(h), order = [h, setdiff(order,h)]; end
    end
    ref = [];
    for i = order
        S = load(fullfile(files(i).folder, files(i).name));
        if ~all(isfield(S,{'Y_pacmap','labels','time_window_idx','time_window_center'})), continue; end
        Y = double(S.Y_pacmap(:,dimsToUse)); lab = double(S.labels(:));
        twi = double(S.time_window_idx(:)); twc = double(S.time_window_center(:));
        sv = unique(lab); sv = sv(isfinite(sv)); if numel(sv)~=2, continue; end
        nogoVal = sv(sv~=goVal); nogoVal = nogoVal(1);
        Cs = makeLandmarks(Y, lab, twi, twc, goVal, nogoVal);
        if sum(all(isfinite(Cs),2)) >= 10
            ref.Cstack = Cs; ref.name = files(i).name; return;
        end
    end
    error('No usable reference session for alignment.');
end

function Cstack = makeLandmarks(Y3, labels, twi, twc, goVal, nogoVal)
    uW = sort(unique(twi(isfinite(twi)))); nW = numel(uW);
    trajGo = nan(nW,3); trajNo = nan(nW,3);
    for wi = 1:nW
        w = uW(wi);
        g = (labels==goVal)   & (twi==w);
        n = (labels==nogoVal) & (twi==w);
        if any(g), trajGo(wi,:) = mean(Y3(g,:),1,'omitnan'); end
        if any(n), trajNo(wi,:) = mean(Y3(n,:),1,'omitnan'); end
    end
    Cstack = [trajGo; trajNo];
end
