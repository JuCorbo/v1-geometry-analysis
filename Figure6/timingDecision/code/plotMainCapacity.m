%% MAIN figure: example licking (left) + capacity divergence timing vs behavior by day (right)
clear; clc;
codeDir = fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
anaDir = fileparts(codeDir);
outDir = codeDir;
FS=6; FN='Arial'; LW=1;

% ---- example session lick data ----
dayStr='D3'; sessIdx=3; goVal=45; noGoVal=75;
stimOnsetSec=1.5; plotWinSec=[-1 4]; doSmooth=true; smoothWinSec=0.15;
crDecisionWin=[0 2];
colGo=[0.22 0.37 0.65]; colNo=[0.85 0.25 0.25]; colStim=[0.4 0.4 0.4]; colLast=[0 0 0]; colDiv=[0.10 0.55 0.25];
cacheFile=fullfile(anaDir,'cache',sprintf('exampleLickCache_%s_s%d.mat',dayStr,sessIdx));
assert(exist(cacheFile,'file')==2, 'Example lick cache not found: %s', cacheFile);
S=load(cacheFile); traces=S.traces; lickDt=S.lickDt; thr=S.thr;
isGo=S.isGo; isNoGo=S.isNoGo; isHit=S.isHit; isFA=S.isFA; isCR=S.isCR; expID=S.expID;
fprintf('Loaded example lick data from cache (%s).\n',cacheFile);
[tVec,goOn]=tracesToLickOn(traces(isGo),lickDt,stimOnsetSec,plotWinSec,doSmooth,smoothWinSec,thr);
[~,noOn]  =tracesToLickOn(traces(isNoGo),lickDt,stimOnsetSec,plotWinSec,doSmooth,smoothWinSec,thr);
pGo=mean(goOn,1,'omitnan'); pNo=mean(noOn,1,'omitnan');
semGo=std(goOn,0,1,'omitnan')./sqrt(sum(~isnan(goOn),1)); semNo=std(noOn,0,1,'omitnan')./sqrt(sum(~isnan(noOn),1));
rng(1); [~,~,dlt,dltLo]=bootstrapDeltaCI(goOn,noOn,300,0.05); tUB=findFirstSustained(tVec,dlt,dltLo,0.10,0.20,[0 6]);
crIdxAll=find(isCR); faIdxAll=find(isFA); hitIdxAll=find(isHit); lastCR=nan(numel(crIdxAll),1);
for i=1:numel(crIdxAll), [~,tl]=firstLastLickTimes(traces(crIdxAll(i)),lickDt,stimOnsetSec,crDecisionWin,doSmooth,smoothWinSec,thr); if ~isempty(tl), lastCR(i)=tl(end); end; end
[~,ord]=sort(lastCR,'descend','MissingPlacement','last'); crIdx=crIdxAll(ord); lastCRs=lastCR(ord); medLast=median(lastCRs,'omitnan');

% ---- capacity timing + behavior ----
% Capacity onset is already reported at the window centre by groupSeparationTest
% (onsetOffsetSec set in capacityBootstrap.m); no local shift here.
R=load(fullfile(anaDir,'data','onsetCI.mat')); capOn=R.results.capacity.onsetPoint; capCI=R.results.capacity.ci;
Tb=readtable(fullfile(anaDir,'data','DecisionUB_ALL_BatuhanAngDist.csv'));
days={'D1','D2','D3','D4','D5','D6'}; colGN=[0.90 0.50 0.10]; colCR=[0.55 0.27 0.07]; colCap=[0.20 0.40 0.80];

% ---- FIGURE ----
fig=figure('Color','w','Units','points','Position',[50 50 620 320]);
axA=axes(fig,'Units','normalized','Position',[0.065 0.70 0.32 0.25]); hold(axA,'on');
axB=axes(fig,'Units','normalized','Position',[0.065 0.09 0.32 0.55]); hold(axB,'on');
axT=axes(fig,'Units','normalized','Position',[0.50 0.11 0.46 0.84]); hold(axT,'on');

% avg
patch(axA,[0 3 3 0],[0 0 1 1],[0.96 0.96 0.96],'EdgeColor','none','HandleVisibility','off');
fill(axA,[tVec fliplr(tVec)],[pGo+semGo fliplr(pGo-semGo)],colGo,'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
fill(axA,[tVec fliplr(tVec)],[pNo+semNo fliplr(pNo-semNo)],colNo,'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
pGh=plot(axA,tVec,pGo,'-','Color',colGo,'LineWidth',LW); pNh=plot(axA,tVec,pNo,'-','Color',colNo,'LineWidth',LW);
xline(axA,0,'--','Color',colStim,'LineWidth',0.75,'HandleVisibility','off');
if isfinite(tUB), xline(axA,tUB,'-','Color',colDiv,'LineWidth',0.75,'HandleVisibility','off'); end
if isfinite(medLast), xline(axA,medLast,'-','Color',colLast,'LineWidth',0.75,'HandleVisibility','off'); end
xlim(axA,plotWinSec); ylim(axA,[0 1]); ylabel(axA,'P(lick)','FontName',FN,'FontSize',FS);
hD=plot(axA,NaN,NaN,'-','Color',colDiv,'LineWidth',LW); hC=plot(axA,NaN,NaN,'-','Color',colLast,'LineWidth',LW);
lg=legend(axA,[pGh pNh hD hC],{'Go','NoGo','Go/NoGo divergence','median last CR lick'},'Location','southoutside','NumColumns',2); set(lg,'Box','off','FontName',FN,'FontSize',FS);
set(axA,'FontName',FN,'FontSize',FS,'LineWidth',LW,'TickDir','out','XTickLabel',[]); box(axA,'off');
title(axA,sprintf('Example  %s  %s',dayStr,expID),'FontName',FN,'FontSize',FS,'FontWeight','normal','Interpreter','none');

% raster
patch(axB,[0 3 3 0],[0.5 0.5 1e4 1e4],[0.96 0.96 0.96],'EdgeColor','none','HandleVisibility','off');
y=0;
for i=1:numel(crIdx), y=y+1; drawLickRow(axB,traces(crIdx(i)),lickDt,stimOnsetSec,plotWinSec,doSmooth,smoothWinSec,thr,y,colNo,0.7);
    if isfinite(lastCRs(i)), plot(axB,lastCRs(i),y,'o','MarkerFaceColor',colLast,'MarkerEdgeColor','none','MarkerSize',2,'HandleVisibility','off'); end; end
for i=1:numel(faIdxAll), y=y+1; drawLickRow(axB,traces(faIdxAll(i)),lickDt,stimOnsetSec,plotWinSec,doSmooth,smoothWinSec,thr,y,[0.95 0.55 0.55],0.7); end
y=y+3;
for i=1:numel(hitIdxAll), y=y+1; drawLickRow(axB,traces(hitIdxAll(i)),lickDt,stimOnsetSec,plotWinSec,doSmooth,smoothWinSec,thr,y,colGo,0.7); end
xline(axB,0,'--','Color',colStim,'LineWidth',0.75,'HandleVisibility','off');
if isfinite(tUB), xline(axB,tUB,'-','Color',colDiv,'LineWidth',0.75,'HandleVisibility','off'); end
if isfinite(medLast), xline(axB,medLast,'-','Color',colLast,'LineWidth',0.75,'HandleVisibility','off'); end
xlim(axB,plotWinSec); ylim(axB,[0.5 y+1]);
xlabel(axB,'Time from stimulus onset (s)','FontName',FN,'FontSize',FS); ylabel(axB,'Trial (CR | FA | Go hits)','FontName',FN,'FontSize',FS);
set(axB,'FontName',FN,'FontSize',FS,'LineWidth',LW,'TickDir','out','YTick',[]); box(axB,'off');

% capacity onset: point+CI (left) + horizontal line across days; behavior by day
YLt=[-0.35 3.4]; xCap=1; xSep=2; xDay=3:8;
grp=Tb.BehavCond;
okG=isfinite(Tb.tUpperBound); [pGN,~]=kruskalwallis(Tb.tUpperBound(okG),grp(okG),'off');
okC=isfinite(Tb.medCR_last);   [pCR,~]=kruskalwallis(Tb.medCR_last(okC),grp(okC),'off');
patch(axT,[2.5 8.5 8.5 2.5],[0 0 3 3],[0.96 0.96 0.96],'EdgeColor','none','HandleVisibility','off');
yline(axT,0,'--','Color',colStim,'LineWidth',0.75,'HandleVisibility','off'); yline(axT,3,'--','Color',colStim,'LineWidth',0.75,'HandleVisibility','off');
plot(axT,[xSep xSep],YLt,'-','Color',[0.85 0.85 0.85],'LineWidth',0.75,'HandleVisibility','off');
plot(axT,[xCap xCap],capCI,'-','Color',colCap,'LineWidth',LW,'HandleVisibility','off');
plot(axT,[xCap-0.12 xCap+0.12],[capCI(1) capCI(1)],'-','Color',colCap,'LineWidth',LW,'HandleVisibility','off');
plot(axT,[xCap-0.12 xCap+0.12],[capCI(2) capCI(2)],'-','Color',colCap,'LineWidth',LW,'HandleVisibility','off');
plot(axT,xCap,capOn,'o','MarkerFaceColor',colCap,'MarkerEdgeColor','k','MarkerSize',6,'LineWidth',0.6,'HandleVisibility','off');
plot(axT,[2.5 8.5],[capOn capOn],'--','Color',colCap,'LineWidth',1,'HandleVisibility','off');
text(axT,8.6,capOn,sprintf('capacity\n%.2f [%.2f, %.2f] s',capOn,capCI(1),capCI(2)),'Color',colCap*0.85,'FontName',FN,'FontSize',FS,'VerticalAlignment','middle');
for d=1:6
    vGN=Tb.tUpperBound(strcmp(Tb.BehavCond,days{d})); vGN=vGN(isfinite(vGN));
    vCR=Tb.medCR_last(strcmp(Tb.BehavCond,days{d}));   vCR=vCR(isfinite(vCR));
    drawVBox(axT,xDay(d)-0.20,vGN,colGN,0.30); drawVBox(axT,xDay(d)+0.20,vCR,colCR,0.30);
end
hGN=plot(axT,NaN,NaN,'s','MarkerFaceColor',colGN,'MarkerEdgeColor',colGN,'MarkerSize',5);
hCRb=plot(axT,NaN,NaN,'s','MarkerFaceColor',colCR,'MarkerEdgeColor',colCR,'MarkerSize',5);
lgT=legend(axT,[hGN hCRb],{'Go/NoGo lick divergence','Last CR lick'},'Location','northwest'); set(lgT,'Box','off','FontName',FN,'FontSize',FS);
text(axT,5.5,YLt(2),sprintf('Day effect (Kruskal-Wallis): Go/NoGo p=%.1e; last CR p=%.2f (n.s.)',pGN,pCR),'FontName',FN,'FontSize',FS,'HorizontalAlignment','center','VerticalAlignment','top');
set(axT,'XTick',[xCap xDay],'XTickLabel',[{'capacity'} days],'FontName',FN,'FontSize',FS,'LineWidth',LW,'TickDir','out'); box(axT,'off');
xlim(axT,[0.3 9.6]); ylim(axT,YLt);
ylabel(axT,'Time from stimulus onset (s)','FontName',FN,'FontSize',FS);
xlabel(axT,'Geometric onset   |   Day (Go/NoGo angle narrows D1\rightarrowD6)','FontName',FN,'FontSize',FS);

exportgraphics(fig,fullfile(anaDir,'panels','main','mainFigure_capacity.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(anaDir,'panels','main','mainFigure_capacity.png'),'Resolution',300);
fprintf('Saved mainFigure_capacity  (tUB=%.2f)\n',tUB);

%% helpers
function drawVBox(ax,x,d,col,w)
    d=d(isfinite(d)); if isempty(d), return; end
    q=prctile(d,[25 50 75]); iqr=q(3)-q(1); lo=max(min(d),q(1)-1.5*iqr); hi=min(max(d),q(3)+1.5*iqr);
    scatter(ax,x+(rand(size(d))-0.5)*w*0.8,d,5,col,'filled','MarkerFaceAlpha',0.40,'HandleVisibility','off');
    plot(ax,[x x],[lo hi],'-','Color',col*0.7,'LineWidth',0.75,'HandleVisibility','off');
    rectangle(ax,'Position',[x-w/2 q(1) w q(3)-q(1)],'FaceColor',[col 0.20],'EdgeColor',col,'LineWidth',0.75);
    plot(ax,[x-w/2 x+w/2],[q(2) q(2)],'-','Color',col*0.55,'LineWidth',1.25,'HandleVisibility','off');
end
function drawLickRow(ax,tr,dt,stimOnsetSec,win,doSmooth,smoothWinSec,thr,y,col,h)
    if iscell(tr), tr=tr{1}; end; if isempty(tr), return; end; v=double(tr(:)');
    if doSmooth&&smoothWinSec>0, v=movmean(v,max(1,round(smoothWinSec/dt)),'Endpoints','shrink'); end
    t=(0:numel(v)-1)*dt - stimOnsetSec; idx=t>=win(1)&t<=win(2); t=t(idx); on=v(idx)>thr;
    d=diff([0 on 0]); s=find(d==1); e=find(d==-1)-1;
    for c=1:numel(s), plot(ax,[t(s(c)) t(e(c))],[y y],'-','Color',col,'LineWidth',h,'HandleVisibility','off'); end
end
function [tVec,lickOn]=tracesToLickOn(traces,dt,stimOnsetSec,win,doSmooth,smoothWinSec,thr)
    tVec=[];lickOn=[]; if isempty(traces)||~(isfinite(dt)&&dt>0), return; end
    lens=cellfun(@numel,traces); Lmax=max(lens); if Lmax<2, return; end; M=NaN(numel(traces),Lmax);
    for k=1:numel(traces), v=traces{k}; if isempty(v), continue; end; v=double(v(:)');
        if doSmooth&&smoothWinSec>0, v=movmean(v,max(1,round(smoothWinSec/dt)),'Endpoints','shrink'); end
        M(k,1:numel(v))=v; end
    tAll=(0:Lmax-1)*dt - stimOnsetSec; idx=tAll>=win(1)&tAll<=win(2); if ~any(idx), return; end
    tVec=tAll(idx); lickOn=M(:,idx)>thr;
end
function [tFirst,tLast]=firstLastLickTimes(traces,dt,stimOnsetSec,win,doSmooth,smoothWinSec,thr)
    tFirst=[];tLast=[]; if iscell(traces)&&numel(traces)==1, v=traces{1}; else, v=traces; end
    if isempty(v), return; end; v=double(v(:)');
    if doSmooth&&smoothWinSec>0, v=movmean(v,max(1,round(smoothWinSec/dt)),'Endpoints','shrink'); end
    t=(0:numel(v)-1)*dt - stimOnsetSec; idx=t>=win(1)&t<=win(2); if ~any(idx), return; end
    tW=t(idx); on=find(v(idx)>thr); if isempty(on), return; end; tFirst=tW(on(1)); tLast=tW(on(end));
end
function [pGo,pNo,delta,deltaLo]=bootstrapDeltaCI(goOn,noOn,nBoot,alpha)
    tLen=size(goOn,2); pGo=mean(goOn,1,'omitnan'); pNo=mean(noOn,1,'omitnan'); delta=pGo-pNo;
    bd=NaN(nBoot,tLen); nG=size(goOn,1); nN=size(noOn,1);
    for b=1:nBoot, bd(b,:)=mean(goOn(randi(nG,nG,1),:),1,'omitnan')-mean(noOn(randi(nN,nN,1),:),1,'omitnan'); end
    deltaLo=prctile(bd,100*(alpha/2),1);
end
function tStar=findFirstSustained(tVec,delta,deltaLo,deltaMin,holdSec,searchWinSec)
    tStar=NaN; if isempty(tVec), return; end
    idxSearch=tVec>=searchWinSec(1)&tVec<=searchWinSec(2); if ~any(idxSearch), return; end
    ok=deltaLo>0; if isfinite(deltaMin)&&deltaMin>0, ok=ok&(delta>deltaMin); end
    ok=ok&idxSearch; dt=median(diff(tVec)); nHold=max(1,round(holdSec/dt));
    okv=double(ok(:))'; dd=diff([0 okv 0]); s=find(dd==1); e=find(dd==-1)-1; g=find((e-s+1)>=nHold,1,'first');
    if ~isempty(g), tStar=tVec(s(g)); end
end
