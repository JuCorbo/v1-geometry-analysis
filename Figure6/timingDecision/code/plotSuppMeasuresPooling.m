%% SUPPLEMENTARY: (top) other geometric measures vs behavior by day;
%  (bottom) pooling justification - cross-day correlation of capacity's
%  trained-naive difference timecourse.
clear; clc;
codeDir = fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
anaDir = fileparts(codeDir);
outDir = codeDir;
FS=6; FN='Arial'; LW=1;
R=load(fullfile(anaDir,'data','onsetCI.mat')); res=R.results;
% Neural onsets are already reported at the window centre by groupSeparationTest
% (onsetOffsetSec set in bootstrapOnsets.m / capacityBootstrap.m); no local
% shift here.
Tb=readtable(fullfile(anaDir,'data','DecisionUB_ALL_BatuhanAngDist.csv'));
days={'D1','D2','D3','D4','D5','D6'};
colStim=[0.4 0.4 0.4]; colGN=[0.90 0.50 0.10]; colCR=[0.55 0.27 0.07];
meas={'dimDPR','dim',[0.80 0.15 0.30],true; 'radiusMix','radius',[0.10 0.65 0.55],true; 'centDist','centroid',[0.60 0.25 0.65],false};
xMetric=1:3; xDay=5:10; xSep=4; YLt=[-0.35 3.7];

% ---- capacity per-day timecourses for pooling justification ----
tcDir=fullfile(anaDir,'cache','capacityTimecourse');
fN=dir(fullfile(tcDir,'*_BatuhanADNaive_capacity_timecourse.mat'));
fT=dir(fullfile(tcDir,'*_BatuhanAngDist_capacity_timecourse.mat'));
capN=cell(1,6); capT=cell(1,6); tvC=[];
for F=[fN;fT]'
    L=load(fullfile(F.folder,F.name),'capacity_rolling','meta','tvec'); if isempty(tvC), tvC=L.tvec(:)'; end
    d=find(strcmp(days,L.meta.BehavCond),1); if isempty(d), continue; end
    if contains(L.meta.DB,'Naive'), capN{d}=[capN{d};L.capacity_rolling(:)']; else, capT{d}=[capT{d};L.capacity_rolling(:)']; end
end
if isempty(tvC), tvC=res.tvec; end
diffTC=nan(6,numel(tvC));
for d=1:6, diffTC(d,:)=mean(capT{d},1,'omitnan')-mean(capN{d},1,'omitnan'); end
% Sustained window where the divergence lives (naive decays / trained maintained)
postW=[1 3]; wp=tvC>=postW(1)&tvC<=postW(2); tp=tvC(wp);
fitLine=nan(6,sum(wp)); for d=1:6, pf=polyfit(tp,diffTC(d,wp),1); fitLine(d,:)=polyval(pf,tp); end
% Pairwise correlation (R) between day difference-curves
Cmat=corr(diffTC(:,wp)','rows','pairwise'); meanR=mean(Cmat(~eye(6)));

% ================= FIGURE =================
fig=figure('Color','w','Units','points','Position',[50 50 640 470]);
axTop=axes(fig,'Units','normalized','Position',[0.08 0.60 0.88 0.33]); hold(axTop,'on');
axTC =axes(fig,'Units','normalized','Position',[0.08 0.11 0.42 0.34]); hold(axTC,'on');
axCr =axes(fig,'Units','normalized','Position',[0.62 0.11 0.26 0.34]); hold(axCr,'on');

% ---- TOP: other measures + behavior by day ----
patch(axTop,[0.3 10.7 10.7 0.3],[0 0 3 3],[0.96 0.96 0.96],'EdgeColor','none','HandleVisibility','off');
plot(axTop,[0.3 10.7],[0 0],'--','Color',colStim,'LineWidth',0.75,'HandleVisibility','off');
plot(axTop,[0.3 10.7],[3 3],'--','Color',colStim,'LineWidth',0.75,'HandleVisibility','off');
plot(axTop,[xSep xSep],YLt,'-','Color',[0.85 0.85 0.85],'LineWidth',0.75,'HandleVisibility','off');
refs={res.dimDPR.onsetPoint,'dim/rad',[0.5 0.2 0.4]; res.centDist.onsetPoint,'centroid',meas{3,3}};
for r=1:size(refs,1)
    plot(axTop,[xSep+0.3 10.7],[refs{r,1} refs{r,1}],'--','Color',refs{r,3},'LineWidth',0.75,'HandleVisibility','off');
    text(axTop,10.8,refs{r,1},refs{r,2},'Color',refs{r,3}*0.8,'FontName',FN,'FontSize',FS,'VerticalAlignment','middle');
end
for m=1:3
    nm=meas{m,1}; c=meas{m,3}; x=xMetric(m); on=res.(nm).onsetPoint; ci=res.(nm).ci;
    if diff(ci)>1e-3
        plot(axTop,[x x],ci,'-','Color',c,'LineWidth',LW,'HandleVisibility','off');
        plot(axTop,[x-0.12 x+0.12],[ci(1) ci(1)],'-','Color',c,'LineWidth',LW,'HandleVisibility','off');
        plot(axTop,[x-0.12 x+0.12],[ci(2) ci(2)],'-','Color',c,'LineWidth',LW,'HandleVisibility','off');
        plot(axTop,x,on,'o','MarkerFaceColor',c,'MarkerEdgeColor','k','MarkerSize',5,'LineWidth',0.5,'HandleVisibility','off');
    else
        plot(axTop,x,on,'v','MarkerFaceColor',c,'MarkerEdgeColor','k','MarkerSize',6,'LineWidth',0.5,'HandleVisibility','off');
    end
end
for d=1:6
    vGN=Tb.tUpperBound(strcmp(Tb.BehavCond,days{d})); vGN=vGN(isfinite(vGN));
    vCR=Tb.medCR_last(strcmp(Tb.BehavCond,days{d}));   vCR=vCR(isfinite(vCR));
    drawVBox(axTop,xDay(d)-0.20,vGN,colGN,0.30); drawVBox(axTop,xDay(d)+0.20,vCR,colCR,0.30);
end
hGN=plot(axTop,NaN,NaN,'s','MarkerFaceColor',colGN,'MarkerEdgeColor',colGN,'MarkerSize',5);
hCR=plot(axTop,NaN,NaN,'s','MarkerFaceColor',colCR,'MarkerEdgeColor',colCR,'MarkerSize',5);
lg=legend(axTop,[hGN hCR],{'Go/NoGo lick divergence','Last CR lick'},'Location','northwest'); set(lg,'Box','off','FontName',FN,'FontSize',FS);
set(axTop,'XTick',[xMetric xDay],'XTickLabel',[meas(:,2)' days],'FontName',FN,'FontSize',FS,'LineWidth',LW,'TickDir','out'); box(axTop,'off');
xlim(axTop,[0.3 11.6]); ylim(axTop,YLt); ylabel(axTop,'Time from stim onset (s)','FontName',FN,'FontSize',FS);
title(axTop,'Other geometric measures (pooled onset +CI) vs behavior by day','FontName',FN,'FontSize',FS,'FontWeight','normal');
text(axTop,0.13,-0.13,'Geometric','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontName',FN,'FontSize',FS);
text(axTop,0.66,-0.13,'Day','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontName',FN,'FontSize',FS);

% ---- BOTTOM-LEFT: capacity diff timecourses per day ----
cmap=[0.10 0.20 0.60; 0.20 0.50 0.85; 0.30 0.70 0.55; 0.75 0.65 0.15; 0.85 0.40 0.10; 0.65 0.15 0.20];
ymax=max(abs(diffTC(:)))*1.1;
patch(axTC,[postW(1) postW(2) postW(2) postW(1)],[-100 -100 100 100],[0.92 0.95 0.99],'EdgeColor','none','HandleVisibility','off');
for d=1:6, plot(axTC,tvC,diffTC(d,:),'-','Color',cmap(d,:),'LineWidth',LW,'DisplayName',days{d}); end
for d=1:6, plot(axTC,tp,fitLine(d,:),'-','Color',cmap(d,:),'LineWidth',3,'HandleVisibility','off'); end   % thick slope fits
xline(axTC,0,'--','Color',colStim,'LineWidth',0.75,'HandleVisibility','off'); yline(axTC,0,':','Color',[0.6 0.6 0.6],'HandleVisibility','off');
xlim(axTC,[min(tvC) max(tvC)]); ylim(axTC,[-0.2*ymax ymax]);
set(axTC,'FontName',FN,'FontSize',FS,'LineWidth',LW,'TickDir','out'); box(axTC,'off');
xlabel(axTC,'Time from stim onset (s)','FontName',FN,'FontSize',FS); ylabel(axTC,'capacity: trained - naive','FontName',FN,'FontSize',FS);
lg2=legend(axTC,'Location','northeast','NumColumns',2); set(lg2,'Box','off','FontName',FN,'FontSize',FS);
title(axTC,'Per-day difference timecourse (bold = 1-3s fit)','FontName',FN,'FontSize',FS,'FontWeight','normal');

% ---- BOTTOM-MIDDLE: pairwise correlation matrix (R) ----
imagesc(axCr,Cmat,[0 1]); axis(axCr,'ij','square','tight'); colormap(axCr,flipud(gray));
set(axCr,'XTick',1:6,'XTickLabel',days,'YTick',1:6,'YTickLabel',days,'FontName',FN,'FontSize',FS,'LineWidth',LW,'TickDir','out');
for i=1:6, for j=1:6, text(axCr,j,i,sprintf('%.2f',Cmat(i,j)),'HorizontalAlignment','center','FontName',FN,'FontSize',FS-1,'Color',[0.85 0.45 0.1]); end; end
title(axCr,sprintf('Cross-day correlation R (1-3s, mean %.2f)',meanR),'FontName',FN,'FontSize',FS,'FontWeight','normal');

exportgraphics(fig,fullfile(anaDir,'panels','supp','suppFigure_measures_pooling.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(anaDir,'panels','supp','suppFigure_measures_pooling.png'),'Resolution',300);
fprintf('Saved suppFigure_measures_pooling  (R mean %.2f)\n',meanR);

function drawVBox(ax,x,d,col,w)
    d=d(isfinite(d)); if isempty(d), return; end
    q=prctile(d,[25 50 75]); iqr=q(3)-q(1); lo=max(min(d),q(1)-1.5*iqr); hi=min(max(d),q(3)+1.5*iqr);
    scatter(ax,x+(rand(size(d))-0.5)*w*0.8,d,5,col,'filled','MarkerFaceAlpha',0.40,'HandleVisibility','off');
    plot(ax,[x x],[lo hi],'-','Color',col*0.7,'LineWidth',0.75,'HandleVisibility','off');
    rectangle(ax,'Position',[x-w/2 q(1) w q(3)-q(1)],'FaceColor',[col 0.20],'EdgeColor',col,'LineWidth',0.75);
    plot(ax,[x-w/2 x+w/2],[q(2) q(2)],'-','Color',col*0.55,'LineWidth',1.25,'HandleVisibility','off');
end
