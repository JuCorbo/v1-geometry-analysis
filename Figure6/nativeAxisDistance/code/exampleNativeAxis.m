%% exampleNativeAxis.m
%  Figure 6 — panel c, NATIVE-space version. The metric lives in full native
%  neural space; for visualization only, the snapshot is projected to its top-3
%  PCs. Go cloud (grey), NoGo cloud colored by the native Go->NoGo axis position
%  (quantile-mapped, Near->Far), the Go->NoGo axis, and a perpendicular
%  illustrative plane. Camera is oriented to the axis so the gradient is visible.
%
%  Reads the shared windowed cache (getWindowed 'global') at 1.5 s.
clear; clc;

%% ---------------- CONFIG ----------------
scheme         = 'global';
snapshotTime   = 1.5;
goVal          = 45;
sessionsToPlot = ["F17_42221_", "F10_42321_"];
keepFrac       = 0.90;
cmapName       = 'parula';

thisDir  = fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir = pwd; end
anaDir   = fileparts(thisDir);
finalRoot= fileparts(fileparts(anaDir));
addpath(fullfile(finalRoot,'_sharedWindowedData','code'));
dataRoot = fullfile(finalRoot,'_sharedWindowedData','data','trained');
panelDir = fullfile(anaDir,'panels','main'); if ~exist(panelDir,'dir'), mkdir(panelDir); end

files = [];
for d = 1:6, files = [files; dir(fullfile(dataRoot,sprintf('D%d',d),'*_win.mat'))]; end %#ok<AGROW>
names = string({files.name});

fig = figure('Color','w','Position',[100 100 1400 640]);
tl  = tiledlayout(1, numel(sessionsToPlot), 'Padding','compact','TileSpacing','compact');

for si = 1:numel(sessionsToPlot)
    hit = find(contains(names, sessionsToPlot(si)),1,'first');
    assert(~isempty(hit),'Session %s not found', sessionsToPlot(si));
    [X, info] = getWindowed(fullfile(files(hit).folder, files(hit).name), scheme);
    tc = info.tvec_centre(:); [~,wk] = min(abs(tc - snapshotTime));
    Xs = squeeze(X(:,:,wk));                    % [nTrials x nCells] native
    L  = double(info.Visual_Stim(:));
    nogoVal = unique(L(L~=goVal)); nogoVal = nogoVal(1);

    % native Go->NoGo axis position (the metric) -> quantile color
    goCn = mean(Xs(L==goVal,:),1); nogoCn = mean(Xs(L==nogoVal,:),1);
    un = (nogoCn-goCn); un = un/norm(un);
    axAll = (Xs - goCn) * un';

    % top-3 PCs for display
    [~, score] = pca(Xs);                        % score is centered
    Y3 = score(:,1:3);
    goC = mean(Y3(L==goVal,:),1); nogoC = mean(Y3(L==nogoVal,:),1);
    u3 = (nogoC-goC); u3 = u3/norm(u3);

    idxGo = (L==goVal); idxNo = (L==nogoVal);
    Yno = Y3(idxNo,:); axNo = axAll(idxNo); qNo = tiedrank(axNo)/numel(axNo);
    if keepFrac<1
        [~,ord]=sort(axNo,'ascend'); nKeep=max(1,floor(keepFrac*numel(axNo)));
        m=false(size(axNo)); m(ord(1:nKeep))=true; YnoP=Yno(m,:); cP=qNo(m);
    else
        YnoP=Yno; cP=qNo;
    end

    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    xlabel(ax,'PC1'); ylabel(ax,'PC2'); zlabel(ax,'PC3');
    scatter3(ax, Y3(idxGo,1),Y3(idxGo,2),Y3(idxGo,3), 16, 0.35*[1 1 1],'filled','MarkerFaceAlpha',0.5);
    scatter3(ax, YnoP(:,1),YnoP(:,2),YnoP(:,3), 30, cP, 'filled','Marker','^','MarkerFaceAlpha',0.9);

    % Go->NoGo axis (line through the two centroids)
    seg = [goC - 0.4*(nogoC-goC); nogoC + 0.6*(nogoC-goC)];
    plot3(ax, seg(:,1),seg(:,2),seg(:,3), 'k-','LineWidth',1.5);
    % perpendicular illustrative plane at the midpoint (normal = u3)
    mid = 0.5*(goC+nogoC);
    e2 = [0 0 1]; if abs(u3(3))>0.9, e2=[0 1 0]; end
    e2 = e2 - (e2*u3')*u3; e2 = e2/norm(e2); e3 = cross(u3,e2);
    sp = 0.6*max(range(Y3,1));
    [gg,hh] = meshgrid(linspace(-sp,sp,12));
    PP = mid + gg(:)*e2 + hh(:)*e3;
    surf(ax, reshape(PP(:,1),size(gg)), reshape(PP(:,2),size(gg)), reshape(PP(:,3),size(gg)), ...
        'FaceAlpha',0.10,'EdgeAlpha',0.08,'FaceColor',[0 0 0]);

    % tight zoom + camera oriented so the axis is roughly horizontal
    xlim(ax,prctile(Y3(:,1),[3 97])); ylim(ax,prctile(Y3(:,2),[3 97])); zlim(ax,prctile(Y3(:,3),[3 97]));
    axis(ax,'vis3d');
    up0=[0 0 1]; if abs(u3(3))>0.9, up0=[0 1 0]; end
    e2c = up0 - (up0*u3')*u3; e2c=e2c/norm(e2c); e3c=cross(u3,e2c);
    ctr=mean(Y3,1); Rc=4*max(range(Y3,1));
    viewOut = (-e3c)*cosd(20) + u3*sind(20); viewOut=viewOut/norm(viewOut);
    set(ax,'CameraTarget',ctr,'CameraUpVector',e2c,'CameraPosition',ctr+Rc*viewOut,'CameraViewAngleMode','auto');

    colormap(ax, feval(cmapName,256));
    cl=[min(cP) max(cP)]; if ~all(isfinite(cl))||diff(cl)<=0, cl=[0 1]; end
    caxis(ax,cl); cb=colorbar(ax); cb.Ticks=cl; cb.TickLabels={'Near','Far'};
    title(ax, sprintf('%s  (native, top-3 PC)  %.2f s', strrep(char(sessionsToPlot(si)),'_','\_'), tc(wk)), 'Interpreter','tex');
end

exportgraphics(fig, fullfile(panelDir,'example_nativeAxis.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig, fullfile(panelDir,'example_nativeAxis.png'),'Resolution',200);
fprintf('Saved native example panel to %s\n', panelDir);
