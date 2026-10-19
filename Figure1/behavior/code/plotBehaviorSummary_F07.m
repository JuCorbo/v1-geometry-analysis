%% plotBehaviorSummary_F07.m — 2-panel behaviour figure (Fig 1c-d), F07 example
%  Panel 1 (Fig 1c): training of example animal F07 (session-resolved;
%           data/fig1_c_training_exampleF07.csv, staged by makeF07TrainingCSV.m).
%  Panel 2: D1-D6 imaging signal-detection (ROC) trajectory, all trained animals,
%           from cache/behaviorSummary_data.mat.
%  Hit%/FA% on the LEFT y-axis, d' on the RIGHT y-axis in Panel 1.
clear; clc; close all;
codeDir=fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
aDir=fileparts(codeDir); cacheDir=fullfile(aDir,'cache'); dataDir=fullfile(aDir,'data'); panelDir=fullfile(aDir,'panels','main');
if ~exist(panelDir,'dir'), mkdir(panelDir); end

cHit=[0.10 0.40 0.80]; cFA=[0.82 0.20 0.20]; cDp=[0 0 0];
clip=@(p)max(min(p,.99),.01); dpr=@(pH,pF)norminv(clip(pH))-norminv(clip(pF));
gngAng=[90 45 30 25 20 15];

%% ---- Panel-2 data from cache (SQL-free) ----
B=load(fullfile(cacheDir,'behaviorSummary_data.mat'));   % pHitAll,pFAAll,dpAll,animals
pHitAll=B.pHitAll; pFAAll=B.pFAAll; nA=size(pHitAll,2);
mHit=mean(pHitAll,2,'omitnan'); mFA=mean(pFAAll,2,'omitnan');
semHit=std(pHitAll,0,2,'omitnan')./sqrt(sum(~isnan(pHitAll),2));
semFA =std(pFAAll ,0,2,'omitnan')./sqrt(sum(~isnan(pFAAll),2));

%% ---- Panel-1 data from F07 CSV ----
E=readtable(fullfile(dataDir,'fig1_c_training_exampleF07.csv'));
xg=E.session; hR=E.hitPct; fR=E.faPct; dR=E.dprime;

%% ============================ FIGURE ============================
fig=figure('Color','w','Units','centimeters','Position',[2 2 18 7.2]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
YL=[0 100]; YR=[-0.5 4.5];

% ---------- Panel 1: training (F07, session-resolved) ----------
ax1=nexttile(tl); yyaxis(ax1,'left'); hold(ax1,'on');
plot(ax1,xg,hR,'-o','Color',cHit,'MarkerFaceColor',cHit,'MarkerSize',4,'LineWidth',2.2,'DisplayName','Hit %');
plot(ax1,xg,fR,'-o','Color',cFA,'MarkerFaceColor',cFA,'MarkerSize',4,'LineWidth',2.2,'DisplayName','FA %');
ylabel(ax1,'rate (%)'); ylim(ax1,YL); ax1.YColor='k';
yyaxis(ax1,'right');
plot(ax1,xg,dR,'-o','Color',cDp,'MarkerFaceColor',cDp,'MarkerSize',4,'LineWidth',2.2,'DisplayName','d''');
yline(ax1,0,':','Color',[.6 .6 .6],'HandleVisibility','off');
yline(ax1,1.5,':','Color',[.6 .6 .6],'HandleVisibility','off');
ylabel(ax1,'d'''); ylim(ax1,YR); ax1.YColor=cDp;
xlabel(ax1,'training session'); xlim(ax1,[0.5 numel(xg)+0.5]);
title(ax1,'Training — example animal (F07)');
legend(ax1,'Location','south','Orientation','horizontal','Box','off'); grid(ax1,'on');

% ---------- Panel 2: signal-detection (ROC) trajectory, all animals ----------
ax2=nexttile(tl); hold(ax2,'on'); axis(ax2,'square');
gv=linspace(1e-3,1-1e-3,300); [FAg,HITg]=meshgrid(gv,gv); Dfield=norminv(HITg)-norminv(FAg);
contourf(ax2,FAg,HITg,Dfield,0:0.2:4,'LineStyle','none');
gg=linspace(0,1,256)'; colormap(ax2,[1-0.42*gg,1-0.28*gg,1-0.08*gg]); clim(ax2,[0 3.5]);
[C,h]=contour(ax2,FAg,HITg,Dfield,0:1:3,'LineColor',[0.55 0.55 0.55],'LineWidth',0.75);
clabel(C,h,'FontSize',7,'Color',[0.4 0.4 0.4],'LabelSpacing',350);
% paper "best random d'" ceiling (Corbo 2025)
nCeil=128; pGoCeil=mean(pHitAll(:),'omitnan');
CEIL=norminv(clip(binoinv(0.775,nCeil,pGoCeil)/nCeil))-norminv(clip(binoinv(0.225,nCeil,0.5)/nCeil));
faC=linspace(1e-3,1-1e-3,300);
plot(ax2,faC,normcdf(norminv(faC)+CEIL),'--','Color',[0 0 0],'LineWidth',2.2);
text(ax2,0.35,normcdf(norminv(0.35)+CEIL),sprintf('best random d''=%.2f',CEIL), ...
    'Color','k','FontSize',8.5,'FontWeight','bold','BackgroundColor','w','Margin',0.6, ...
    'HorizontalAlignment','center','VerticalAlignment','bottom');
dayCol=turbo(6);
for a=1:nA, plot(ax2,pFAAll(:,a),pHitAll(:,a),'-','Color',[0.45 0.45 0.45 0.30],'LineWidth',1.2); end
for d=1:6, scatter(ax2,pFAAll(d,:),pHitAll(d,:),4,dayCol(d,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.85); end
plot(ax2,mFA,mHit,'-','Color',[0 0 0 0.8],'LineWidth',2);
for d=1:6
    plot(ax2,[mFA(d)-semFA(d) mFA(d)+semFA(d)],[mHit(d) mHit(d)],'-k','LineWidth',1.3);
    plot(ax2,[mFA(d) mFA(d)],[mHit(d)-semHit(d) mHit(d)+semHit(d)],'-k','LineWidth',1.3);
end
for d=1:6
    plot(ax2,mFA(d),mHit(d),'o','MarkerFaceColor',dayCol(d,:),'MarkerEdgeColor','k','MarkerSize',5,'LineWidth',0.75);
    text(ax2,mFA(d)+0.03,mHit(d)-0.03,sprintf('D%d',d),'FontWeight','bold','FontSize',9);
end
xlim(ax2,[0 1]); ylim(ax2,[0 1]); box(ax2,'on');
xlabel(ax2,'false-alarm rate'); ylabel(ax2,'hit rate');
title(ax2,sprintf('Imaging D1–D6 — signal-detection trajectory (n=%d)',nA));
cb=colorbar(ax2); cb.Label.String='sensitivity  d''';

title(tl,'Trained v1 mice — behavioral performance','FontWeight','bold');
for a=[ax1 ax2], set(a,'FontName','Arial','FontSize',6,'LineWidth',1,'Box','off'); end
tl.Title.FontName='Arial'; tl.Title.FontSize=6;
set(findall(fig,'Type','text'),'FontName','Arial','FontSize',6);
set(findall(fig,'Type','colorbar'),'FontName','Arial','FontSize',6);
set(findall(fig,'Type','legend'),'FontName','Arial','FontSize',6);

base=fullfile(panelDir,'fig1_cd_behaviorSummary');
savefig(fig,[base '.fig']); print(fig,[base '.png'],'-dpng','-r300');
try, exportgraphics(fig,[base '.pdf'],'ContentType','vector'); catch, end
fprintf('DONE -> %s.(png/pdf/fig)\n',base);
