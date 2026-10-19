%% plotLickerTypes.m  — Figure 6 supplement: detection vs decision lickers
%  Two example Go/NoGo lick-probability timecourses (a detection licker, whose Go
%  licking rises right at stimulus onset, and a decision licker, whose Go licking
%  rises late, timed to the reward in the third second) + the population
%  distribution of the strategy timing go_t50 (time at which the Go lick
%  probability reaches 50% of its peak above baseline).
%
%  Inputs: go_t50 + per-session lick threshold from ../data/DecisionUB CSV;
%  example lick timecourses from the cached full lick traces (timing module).
clear; clc;
codeDir=fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
anaDir=fileparts(codeDir);
csv=fullfile(anaDir,'data','DecisionUB_ALL_BatuhanAngDist.csv');
lickCacheRoot=fullfile(fileparts(fileparts(anaDir)),'Figure5','behavioralControls','cache');
lickDirs={fullfile(lickCacheRoot,'lickTracesFull'), fullfile(lickCacheRoot,'lickTracesFull_D2')};
suppDir=fullfile(anaDir,'panels','supp'); if ~exist(suppDir,'dir'), mkdir(suppDir); end
statsDir=fullfile(anaDir,'stats'); if ~exist(statsDir,'dir'), mkdir(statsDir); end

stimOnsetSec=1.5; smoothWinSec=0.15; plotWin=[-1 4]; rewardWin=[2 3];  % 3rd second of the 3 s stimulus
goVal=45; colGo=[0.15 0.35 0.85]; colNo=[0.85 0.2 0.2];
colDet=[0.20 0.60 0.55]; colDec=[0.90 0.55 0.20];   % detection / decision

T=readtable(csv);
t50=T.go_t50; amp=T.go_amp; okAll=isfinite(t50);

% file map exp -> cached lick traces
fmap=containers.Map;
for d=1:numel(lickDirs)
  if ~isfolder(lickDirs{d}), continue; end
  ff=dir(fullfile(lickDirs{d},'*_BatuhanAngDist_lickTraces.mat'));
  for f=1:numel(ff), id=extractBefore(ff(f).name,'_BatuhanAngDist_lickTraces.mat');
    if ~isempty(id)&&~isKey(fmap,id), fmap(id)=fullfile(ff(f).folder,ff(f).name); end; end
end

% detection example = the session used in the main figure (C05_121220, go_t50~0.74);
% decision example = a clean late-t50 session (auto).
detExp='C05_121220';   % detection = main-figure example (D3, go_t50~0.74)
decExp='F17_42221_';   % decision  = distance-thread example (D3, go_t50~1.86)
cand=find(okAll & amp>0.35 & cellfun(@(e) isKey(fmap,e), T.Experiment));
[~,ord]=sort(t50(cand));
detIdx=find(strcmp(T.Experiment,detExp),1);
if isempty(detIdx), detIdx=cand(ord(max(1,round(0.10*numel(ord))))); end   % fallback
decIdx=find(strcmp(T.Experiment,decExp),1);
if isempty(decIdx), decIdx=cand(ord(min(numel(ord),round(0.92*numel(ord))))); end   % fallback
exIdx=[detIdx decIdx];
exLab={'Detection licker','Decision licker'};

fig=figure('Color','w','Position',[70 90 1400 420]);
tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k=1:2
  r=exIdx(k); exp=T.Experiment{r}; thr=T.lickThresh(r);
  L=load(fmap(exp),'lickTrace','Visual_Stim','Lick_dt');
  vs=double(L.Visual_Stim); dt=median(double(L.Lick_dt(isfinite(double(L.Lick_dt)))));
  [tVec,goOn]=tracesToLickOn(L.lickTrace(vs==goVal),dt,stimOnsetSec,plotWin,true,smoothWinSec,thr);
  [~,noOn]   =tracesToLickOn(L.lickTrace(vs~=goVal),dt,stimOnsetSec,plotWin,true,smoothWinSec,thr);
  pGo=mean(goOn,1,'omitnan'); pNo=mean(noOn,1,'omitnan');
  ax=nexttile(tl,k); hold(ax,'on'); box(ax,'off');
  patch(ax,[rewardWin(1) rewardWin(2) rewardWin(2) rewardWin(1)],[0 0 1 1],[0.9 0.9 0.6],'EdgeColor','none','FaceAlpha',0.4);
  xline(ax,0,'--','Color',[0.5 0.5 0.5]);
  plot(ax,tVec,pGo,'-','Color',colGo,'LineWidth',2);
  plot(ax,tVec,pNo,'-','Color',colNo,'LineWidth',2);
  if isfinite(t50(r)), xline(ax,t50(r),'-','Color',colGo,'LineWidth',1.3); end
  ylim(ax,[0 1]); xlim(ax,plotWin); xlabel(ax,'Time from stimulus onset (s)');
  if k==1, ylabel(ax,'P(lick)'); end
  title(ax,sprintf('%s  %s  (go\\_t50=%.2f s)', char(9312+k-1), exLab{k}, t50(r)),'Interpreter','tex');
  if k==1
    text(ax,0.98,0.95,'Go','Color',colGo,'Units','normalized','HorizontalAlignment','right','FontWeight','bold');
    text(ax,0.98,0.85,'NoGo','Color',colNo,'Units','normalized','HorizontalAlignment','right','FontWeight','bold');
    text(ax,mean(rewardWin),0.02,'reward window','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7,'Color',[0.5 0.5 0.3]);
  end
end

% distribution of go_t50 across sessions
ax=nexttile(tl,3); hold(ax,'on'); box(ax,'off');
v=t50(okAll);
histogram(ax,v,20,'FaceColor',0.6*[1 1 1],'EdgeColor','none');
yl=ylim(ax);
plot(ax,[v v]',repmat([0;-0.04*yl(2)],1,numel(v)),'-','Color',[0 0 0 0.25],'LineWidth',0.5);  % rug
for k=1:2, xline(ax,t50(exIdx(k)),'-','Color',colGo,'LineWidth',2); ...
    text(ax,t50(exIdx(k)),yl(2),sprintf(' %s',char(9312+k-1)),'Color',colGo,'VerticalAlignment','top'); end
% split detection vs decision at the antimode (density minimum) of the bimodal go_t50
[fk,xk]=ksdensity(v); reg=xk>=0.9 & xk<=1.4;
if any(reg), xr=xk(reg); [~,mi]=min(fk(reg)); splitThr=xr(mi); else, splitThr=1.0; end
pctDet=100*mean(v<splitThr); pctDec=100*mean(v>=splitThr);
xline(ax,splitThr,':','Color',[0.3 0.3 0.3],'LineWidth',1.2);
xline(ax,median(v),'k--','LineWidth',1);
xlabel(ax,'go\_t50 (s): detection \leftarrow  \rightarrow decision'); ylabel(ax,'sessions');
title(ax,sprintf('%s  Strategy-timing distribution (n=%d, median=%.2f s)',char(9314),numel(v),median(v)),'Interpreter','tex');
ylim(ax,[yl(1) yl(2)*1.5]);                       % headroom so the pie clears the bars
% pie inset (top-right): fraction of sessions in each licker type
drawnow; pos=ax.Position; pw=0.085;
axp=axes('Position',[pos(1)+pos(3)*0.60, pos(2)+pos(4)-pw-0.028, pw, pw]);
ph=pie(axp,[pctDet pctDec], {sprintf('detection %.0f%%',pctDet), sprintf('decision %.0f%%',pctDec)});
set(ph(1),'FaceColor',colDet,'EdgeColor','w'); set(ph(3),'FaceColor',colDec,'EdgeColor','w');
ph(2).Color=colDet; ph(2).FontSize=7; ph(2).FontWeight='bold';
ph(4).Color=colDec; ph(4).FontSize=7; ph(4).FontWeight='bold';

exportgraphics(fig,fullfile(suppDir,'fig6_supp_lickerTypes.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(suppDir,'fig6_supp_lickerTypes.png'),'Resolution',200);

fid=fopen(fullfile(statsDir,'lickerTypes_stats.txt'),'w');
fprintf(fid,'Detection vs decision lickers (Figure 6 supp)\n');
fprintf(fid,'go_t50 = time Go P(lick) reaches 50%% of peak-above-baseline (s from stim onset)\n');
fprintf(fid,'distribution: n=%d, median=%.3f, IQR[%.3f, %.3f], range[%.3f, %.3f]\n', ...
  numel(v),median(v),prctile(v,25),prctile(v,75),min(v),max(v));
fprintf(fid,'examples: detection %s go_t50=%.3f | decision %s go_t50=%.3f\n', ...
  T.Experiment{detIdx},t50(detIdx),T.Experiment{decIdx},t50(decIdx));
fprintf(fid,'split at antimode go_t50=%.3f s -> detection %.0f%% (n=%d) / decision %.0f%% (n=%d)\n', ...
  splitThr, pctDet, sum(v<splitThr), pctDec, sum(v>=splitThr));
fclose(fid);
fprintf('detection=%s (%.2f)  decision=%s (%.2f)  | median go_t50=%.2f s (n=%d)\n', ...
  T.Experiment{detIdx},t50(detIdx),T.Experiment{decIdx},t50(decIdx),median(v),numel(v));
fprintf('Saved fig6_supp_lickerTypes to %s\n',suppDir);

%% helpers (lick binarisation, from recomputeDivergence_fromCache)
function [tVec,lickOn]=tracesToLickOn(traces,dt,stimOnsetSec,win,doSmooth,smoothWinSec,thr)
    tVec=[];lickOn=[]; if isempty(traces)||~(isfinite(dt)&&dt>0), return; end
    lens=cellfun(@numel,traces); Lmax=max(lens); if Lmax<2, return; end
    M=NaN(numel(traces),Lmax);
    for k=1:numel(traces), vv=traces{k}; if isempty(vv), continue; end; vv=double(vv(:)');
        if doSmooth&&smoothWinSec>0, vv=movmean(vv,max(1,round(smoothWinSec/dt)),'Endpoints','shrink'); end
        M(k,1:numel(vv))=vv; end
    tAll=(0:Lmax-1)*dt - stimOnsetSec; idx=tAll>=win(1)&tAll<=win(2); if ~any(idx), return; end
    tVec=tAll(idx); lickOn=M(:,idx)>thr;
end
