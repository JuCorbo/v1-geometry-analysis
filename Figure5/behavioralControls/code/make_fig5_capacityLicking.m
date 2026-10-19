%% make_fig5_capacityLicking.m  —  Figure 5 (behavioural controls: licking tertiles)
%  - capacity (T1/T3 tertiles + randThirdAvg) from the capacityTimecourse caches
%  - lick traces from the cached full export (lickTracesFull; un-truncated 13.5 s)
%  - tertile metric LickAUCz_1_2 from the per-trial AUC csv ([1,2] s window)
%  - session->day from the shared windowed cache
% Panels: A example sessions | B day aggregate (+ Naive 1/3 line) | C scatter | D % boxes.
clear; clc;

%% ---------------- paths ----------------
codeDir = fileparts(mfilename('fullpath'));
anaDir  = fileparts(codeDir);                                   % behavioralControls
figDir  = fileparts(fileparts(anaDir));                         % Figures_final
capDir  = fullfile(figDir,'Figure4','manifoldCapacity','cache','capacityTimecourse');
lickDir = fullfile(anaDir,'cache','lickTracesFull');
aucCsv  = fullfile(anaDir,'cache','LickAUCz_perTrial.csv');
sharedTr= fullfile(figDir,'_sharedWindowedData','data','trained');
sharedNv= fullfile(figDir,'_sharedWindowedData','data','naive');
panelDir= fullfile(anaDir,'panels','main'); if ~exist(panelDir,'dir'), mkdir(panelDir); end
dataDir = fullfile(anaDir,'data');  if ~exist(dataDir,'dir'),  mkdir(dataDir);  end
statsDir= fullfile(anaDir,'stats'); if ~exist(statsDir,'dir'), mkdir(statsDir); end

%% ---------------- config (matches original) ----------------
DB='BatuhanAngDist'; dbNaive='BatuhanADNaive'; lickCol='LickAUCz_1_2';
BehavCond={'D1','D2','D3','D4','D5','D6'};
goVal=45; quantWin=[1 2]; stimOnsetSec=1.5; tPlotRel=[-1 4];
doSmooth=true; smoothWinSec=0.15;
exampleDay='D2'; exampleSessions={'B10_111020','F10_42221_'};
colT1=[0.85 0.4 0.15]; colT3=[0.10 0.85 0.2]; lsGo='-'; lsNo='--';
lwCapSess=1.0; lwCapMean=2.8; lwCapEx=2.6; lwLickMean=2.0; lwLickEx=2.0;
alphaSess=0.22; alphaWin=0.10; alphaDiff=0.18; colWin=[0 0 0]; colDiff=[0 0 0];
scatterDayColors=lines(numel(BehavCond)); scatterCrossLW=2.8;
exDayIdx=find(strcmp(BehavCond,exampleDay),1);

%% ---------------- day -> sessions (from shared cache) ----------------
daySess=cell(6,1); naiveDaySess=cell(6,1);
for d=1:6
  ff=dir(fullfile(sharedTr,sprintf('D%d',d),'*_win.mat'));
  daySess{d}=erase(erase({ff.name},[DB '_']),'_win.mat');
  gg=dir(fullfile(sharedNv,sprintf('D%d',d),'*_win.mat'));
  naiveDaySess{d}=erase(erase({gg.name},[dbNaive '_']),'_win.mat');
end

%% ---------------- AUC map (Experiment|Trial -> LickAUCz_1_2) ----------------
Q=readtable(aucCsv);
aucMap=containers.Map('KeyType','char','ValueType','double');
for i=1:height(Q), aucMap([char(string(Q.Experiment(i))) '|' num2str(Q.Trial(i))])=Q.LickAUCz_1_2(i); end

%% ============================================================
% (A) EXAMPLE SESSIONS
%% ============================================================
figA=figure('Color','w','Name',sprintf('Fig A - %s examples',exampleDay),'Position',[60 60 760 620]);
tlA=tiledlayout(figA,2,numel(exampleSessions),'Padding','compact','TileSpacing','compact');
title(tlA,sprintf('Fig A - Example sessions (%s): capacity (top) and licking tertiles (bottom)',exampleDay),'Interpreter','none');
for k=1:numel(exampleSessions)
  sess=exampleSessions{k};
  if ~hasTierCap(capDir,DB,lickCol,sess), warning('example %s missing tier cap, skip',sess); continue; end
  [capT1,capT3,tvec]=loadTierCap(capDir,DB,lickCol,sess);
  axCap=nexttile(tlA,k); hold(axCap,'on'); title(axCap,sprintf('%s - Capacity',sess),'Interpreter','none');
  addDiffPatch(axCap,tvec,capT1,capT3,quantWin,alphaDiff,colDiff);
  plot(axCap,tvec,capT1,'-','Color',colT1,'LineWidth',lwCapEx);
  plot(axCap,tvec,capT3,'-','Color',colT3,'LineWidth',lwCapEx);
  xline(axCap,0,':k','LineWidth',1); xlabel(axCap,'Time (s rel. stim)'); ylabel(axCap,'Capacity');
  box(axCap,'off'); axis(axCap,'tight'); ylim(axCap,[0 0.15]); addWindowShade(axCap,quantWin,alphaWin,colWin);
  pctLabel(axCap,tvec,capT1,capT3,quantWin,colDiff);
  [tlick,mGoT1,mGoT3,mNoT1,mNoT3]=meanLickTierGoNo(lickDir,aucMap,DB,sess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec);
  axL=nexttile(tlA,numel(exampleSessions)+k); hold(axL,'on'); title(axL,sprintf('%s - Licking',sess),'Interpreter','none');
  addDiffPatch(axL,tlick(:),mGoT1(:),mGoT3(:),quantWin,alphaDiff,colDiff);
  addDiffPatch(axL,tlick(:),mNoT1(:),mNoT3(:),quantWin,alphaDiff,colDiff);
  plot(axL,tlick,mGoT1,lsGo,'Color',colT1,'LineWidth',lwLickEx);
  plot(axL,tlick,mGoT3,lsGo,'Color',colT3,'LineWidth',lwLickEx);
  plot(axL,tlick,mNoT1,lsNo,'Color',colT1,'LineWidth',lwLickEx);
  plot(axL,tlick,mNoT3,lsNo,'Color',colT3,'LineWidth',lwLickEx);
  xline(axL,0,':k','LineWidth',1); xlabel(axL,'Time (s rel. stim)'); ylabel(axL,'Licking (a.u.)');
  box(axL,'off'); axis(axL,'tight'); addWindowShade(axL,quantWin,alphaWin,colWin);
  lickLabel(axL,tlick,mGoT1,mGoT3,mNoT1,mNoT3,quantWin,colDiff);
end
saveFig(figA,fullfile(panelDir,'fig5a_examples'));

%% ============================================================
% (B) DAY AGGREGATE (+ Naive 1/3 line)
%% ============================================================
figB=figure('Color','w','Name',sprintf('Fig B - %s aggregate',exampleDay),'Position',[60 60 720 640]);
tlB=tiledlayout(figB,2,1,'Padding','compact','TileSpacing','compact');
title(tlB,sprintf('Fig B - %s: all sessions (capacity top, licking bottom)',exampleDay),'Interpreter','none');
axCapB=nexttile(tlB,1); hold(axCapB,'on'); axLickB=nexttile(tlB,2); hold(axLickB,'on');
capT1_all=[];capT3_all=[];tvec_ref=[]; lGoT1=[];lGoT3=[];lNoT1=[];lNoT3=[];tlick_ref=[];
sList=daySess{exDayIdx};
for si=1:numel(sList)
  sess=char(sList{si}); if ~hasTierCap(capDir,DB,lickCol,sess), continue; end
  [capT1,capT3,tvec]=loadTierCap(capDir,DB,lickCol,sess);
  if isempty(tvec_ref), tvec_ref=tvec; end
  [~,c1A,c3A]=alignToReference(tvec_ref,tvec,capT1,capT3);
  p=plot(axCapB,tvec_ref,c1A,'-','LineWidth',lwCapSess); p.Color=[colT1 alphaSess];
  p=plot(axCapB,tvec_ref,c3A,'-','LineWidth',lwCapSess); p.Color=[colT3 alphaSess];
  capT1_all=[capT1_all;c1A(:)']; capT3_all=[capT3_all;c3A(:)'];
  [tlick,mGoT1,mGoT3,mNoT1,mNoT3]=meanLickTierGoNo(lickDir,aucMap,DB,sess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec);
  if isempty(tlick), continue; end
  if isempty(tlick_ref), tlick_ref=tlick(:); end
  [~,g1,g3]=alignToReference(tlick_ref,tlick(:),mGoT1(:),mGoT3(:));
  [~,n1,n3]=alignToReference(tlick_ref,tlick(:),mNoT1(:),mNoT3(:));
  lGoT1=[lGoT1;g1(:)']; lGoT3=[lGoT3;g3(:)']; lNoT1=[lNoT1;n1(:)']; lNoT3=[lNoT3;n3(:)'];
end
capT1_m=mean(capT1_all,1,'omitnan')'; capT3_m=mean(capT3_all,1,'omitnan')';
addDiffPatch(axCapB,tvec_ref,capT1_m,capT3_m,quantWin,alphaDiff,colDiff);
pCap(1)=plot(axCapB,tvec_ref,capT1_m,'-','Color',colT1,'LineWidth',lwCapMean);
pCap(2)=plot(axCapB,tvec_ref,capT3_m,'-','Color',colT3,'LineWidth',lwCapMean);
% naive random-thirds line (same day)
capNaive_all=[];
for si=1:numel(naiveDaySess{exDayIdx})
  s=char(naiveDaySess{exDayIdx}{si});
  f=fullfile(capDir,sprintf('%s_%s_randThirdAvg_capacity_timecourse.mat',s,dbNaive));
  if ~isfile(f), continue; end
  Sn=load(f,'capacity_rolling_randThirdAvg','tvec');
  [~,cNA,~]=alignToReference(tvec_ref,Sn.tvec(:),Sn.capacity_rolling_randThirdAvg(:),Sn.capacity_rolling_randThirdAvg(:));
  capNaive_all=[capNaive_all;cNA(:)'];
end
if ~isempty(capNaive_all)
  capNaive_m=mean(capNaive_all,1,'omitnan')';
  pCap(3)=plot(axCapB,tvec_ref,capNaive_m,':','Color',[0 0 0],'LineWidth',lwCapMean);
end
xline(axCapB,0,':k','LineWidth',1); xlabel(axCapB,'Time (s rel. stim)'); ylabel(axCapB,'Capacity');
box(axCapB,'off'); axis(axCapB,'tight'); ylim(axCapB,[0 0.2]); addWindowShade(axCapB,quantWin,alphaWin,colWin);
pctLabel(axCapB,tvec_ref,capT1_m,capT3_m,quantWin,colDiff);
if numel(pCap)>=3, legend(pCap,{'T1 cap (mean)','T3 cap (mean)','Naive (1/3 trials)'},'Location','bestoutside');
else, legend(pCap,{'T1 cap (mean)','T3 cap (mean)'},'Location','bestoutside'); end
lGoT1_m=mean(lGoT1,1,'omitnan')'; lGoT3_m=mean(lGoT3,1,'omitnan')';
lNoT1_m=mean(lNoT1,1,'omitnan')'; lNoT3_m=mean(lNoT3,1,'omitnan')';
addDiffPatch(axLickB,tlick_ref,lGoT1_m,lGoT3_m,quantWin,alphaDiff,colDiff);
addDiffPatch(axLickB,tlick_ref,lNoT1_m,lNoT3_m,quantWin,alphaDiff,colDiff);
pL(1)=plot(axLickB,tlick_ref,lGoT1_m,lsGo,'Color',colT1,'LineWidth',lwLickMean);
pL(2)=plot(axLickB,tlick_ref,lGoT3_m,lsGo,'Color',colT3,'LineWidth',lwLickMean);
pL(3)=plot(axLickB,tlick_ref,lNoT1_m,lsNo,'Color',colT1,'LineWidth',lwLickMean);
pL(4)=plot(axLickB,tlick_ref,lNoT3_m,lsNo,'Color',colT3,'LineWidth',lwLickMean);
xline(axLickB,0,':k','LineWidth',1); xlabel(axLickB,'Time (s rel. stim)'); ylabel(axLickB,'Licking (a.u.)');
box(axLickB,'off'); axis(axLickB,'tight'); addWindowShade(axLickB,quantWin,alphaWin,colWin);
lickLabel(axLickB,tlick_ref,lGoT1_m,lGoT3_m,lNoT1_m,lNoT3_m,quantWin,colDiff);
legend(pL,{'GO T1','GO T3','NOGO T1','NOGO T3'},'Location','bestoutside');
saveFig(figB,fullfile(panelDir,'fig5b_dayAggregate'));

%% ============================================================
% (C) SCATTER: dCapacity vs dLickAll (abs, all days)
%% ============================================================
[dCap,dLick,dDay]=collectScatter(capDir,lickDir,aucMap,DB,lickCol,daySess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec,quantWin);
figC=figure('Color','w','Name','Fig C - dCap vs dLick','Position',[60 60 640 600]); axC=axes(figC); hold(axC,'on');
title(axC,sprintf('Fig C - All days: \\DeltaCapacity vs \\DeltaLickAll (%g-%g s)',quantWin(1),quantWin(2)));
xlabel(axC,sprintf('mean \\DeltaLickAll (T3-T1) in %g-%g s',quantWin(1),quantWin(2)));
ylabel(axC,sprintf('mean \\DeltaCapacity (T3-T1) in %g-%g s',quantWin(1),quantWin(2)));
meanP=gobjects(6,1);
for d=1:6
  idx=(dDay==d)&isfinite(dCap)&isfinite(dLick); if ~any(idx), continue; end
  scatter(axC,dLick(idx),dCap(idx),35,scatterDayColors(d,:),'filled','MarkerFaceAlpha',0.22);
  meanP(d)=plot(axC,mean(dLick(idx),'omitnan'),mean(dCap(idx),'omitnan'),'x','Color',scatterDayColors(d,:),'MarkerSize',12,'LineWidth',scatterCrossLW);
end
ok=isfinite(dCap)&isfinite(dLick); r=NaN; pv=NaN;
if nnz(ok)>=3
  x=dLick(ok); y=dCap(ok); pf=polyfit(x,y,1); xf=linspace(min(x),max(x),100);
  plot(axC,xf,polyval(pf,xf),'k-','LineWidth',2); [r,pv]=corr(x,y,'Type','Pearson');
  text(axC,0.05,0.95,sprintf('Pearson r=%.2f, p=%.3g (n=%d)',r,pv,nnz(ok)),'Units','normalized','VerticalAlignment','top');
end
xline(axC,0,':k'); yline(axC,0,':k'); box(axC,'off'); axis(axC,'square');
legend(axC,meanP,BehavCond,'Location','bestoutside');
saveFig(figC,fullfile(panelDir,'fig5c_scatter'));

%% ============================================================
% (D) % CHANGE BOXES: Lick T3-T1 | Cap T3-T1 | Cap Trained-Naive(thirds)
%% ============================================================
[cap1m,cap3m,lick1m,lick3m,~]=collectWinMeans(capDir,lickDir,aucMap,DB,lickCol,daySess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec,quantWin);
denCap=(abs(cap3m)+abs(cap1m))/2; denLick=(abs(lick3m)+abs(lick1m))/2;
pctCap=100*(cap3m-cap1m)./denCap; pctLick=100*(lick3m-lick1m)./denLick;
pctCap(denCap==0)=NaN; pctLick(denLick==0)=NaN;
okPct=isfinite(pctCap)&isfinite(pctLick); pctCap=pctCap(okPct); pctLick=pctLick(okPct);
% trained-vs-naive random thirds, per day
pctTrNv=[];
for d=1:6
  nm=[];
  for si=1:numel(naiveDaySess{d})
    s=char(naiveDaySess{d}{si}); f=fullfile(capDir,sprintf('%s_%s_randThirdAvg_capacity_timecourse.mat',s,dbNaive));
    if ~isfile(f), continue; end; S=load(f,'capacity_rolling_randThirdAvg','tvec');
    idxW=S.tvec(:)>=quantWin(1)&S.tvec(:)<=quantWin(2); nm(end+1,1)=mean(S.capacity_rolling_randThirdAvg(idxW),'omitnan'); %#ok<AGROW>
  end
  naiveDayMean=mean(nm,'omitnan');
  for si=1:numel(daySess{d})
    s=char(daySess{d}{si}); f=fullfile(capDir,sprintf('%s_%s_randThirdAvg_capacity_timecourse.mat',s,DB));
    if ~isfile(f), continue; end; S=load(f,'capacity_rolling_randThirdAvg','tvec');
    idxW=S.tvec(:)>=quantWin(1)&S.tvec(:)<=quantWin(2); tw=mean(S.capacity_rolling_randThirdAvg(idxW),'omitnan');
    den=(abs(tw)+abs(naiveDayMean))/2; if den==0, continue; end
    pctTrNv(end+1,1)=100*(tw-naiveDayMean)/den; %#ok<AGROW>
  end
end
figD=figure('Color','w','Name','Fig D - percent change boxes','Position',[60 60 560 520]); axD=axes(figD); hold(axD,'on');
m=max([numel(pctLick),numel(pctCap),numel(pctTrNv)]); data=nan(m,3);
data(1:numel(pctLick),1)=pctLick; data(1:numel(pctCap),2)=pctCap; data(1:numel(pctTrNv),3)=pctTrNv;
boxplot(axD,data,'Labels',{'Licking T3-T1','Capacity T3-T1','Capacity Tr-Nv'},'Whisker',1.5);
scatter(axD,1+0.05*randn(size(pctLick)),pctLick,30,colT3,'filled','MarkerFaceAlpha',0.35);
scatter(axD,2+0.05*randn(size(pctCap)),pctCap,30,colDiff,'filled','MarkerFaceAlpha',0.35);
scatter(axD,3+0.05*randn(size(pctTrNv)),pctTrNv,30,[0.2 0.4 0.8],'filled','MarkerFaceAlpha',0.35);
yline(axD,0,'k--'); ylabel(axD,'Symmetric % change'); title(axD,sprintf('Symmetric %% change [%.0f %.0f] s',quantWin(1),quantWin(2)));
box(axD,'off');
saveFig(figD,fullfile(panelDir,'fig5d_pctBoxes'));

%% ============================================================
% stats + plot-data CSV
%% ============================================================
cap_ok=dCap(isfinite(dCap)); lick_ok=dLick(isfinite(dLick));
[pC0,~,sC]=signrank(cap_ok,0); [pL0,~,sL]=signrank(lick_ok,0);
% Hedges g_av on capacity window means (T3 vs T1)
x1=cap1m(okPct); x2=cap3m(okPct); N=numel(x1); Dd=x2-x1; sd_av=sqrt((std(x1,'omitnan')^2+std(x2,'omitnan')^2)/2);
g_av=(1-3/(4*N-9))*(mean(Dd,'omitnan')/sd_av);
fid=fopen(fullfile(statsDir,'fig5_capacityLicking_stats.txt'),'w');
fprintf(fid,'FIGURE 5 main - licking tertiles (capacity vs licking), trained only.\n');
fprintf(fid,'Tertiles on %s ([1,2]s lick AUC); capacity and lick traces from the cached intermediates.\n\n',lickCol);
fprintf(fid,'Panel C scatter (all days, n=%d sessions): Pearson r=%.3f, p=%.4g\n',nnz(ok),r,pv);
fprintf(fid,'  signrank dCapacity vs 0: median=%.4g p=%.3g (z=%.3f n=%d)\n',median(cap_ok,'omitnan'),pC0,sC.zval,numel(cap_ok));
fprintf(fid,'  signrank dLickAll  vs 0: median=%.4g p=%.3g (z=%.3f n=%d)\n\n',median(lick_ok,'omitnan'),pL0,sL.zval,numel(lick_ok));
fprintf(fid,'Panel D symmetric %% change [1,2]s (median):\n');
fprintf(fid,'  Licking T3-T1:  median=%.1f%% (n=%d)\n',median(pctLick,'omitnan'),numel(pctLick));
fprintf(fid,'  Capacity T3-T1: median=%.1f%% (n=%d)  Hedges g_av=%.3f (N=%d)\n',median(pctCap,'omitnan'),numel(pctCap),g_av,N);
fprintf(fid,'  Capacity Tr-Nv: median=%.1f%% (n=%d)\n',median(pctTrNv,'omitnan'),numel(pctTrNv));
fclose(fid); type(fullfile(statsDir,'fig5_capacityLicking_stats.txt'));
% plot-data csv (scatter points)
T=table(dLick(ok),dCap(ok),dDay(ok),'VariableNames',{'dLickAll_T3mT1','dCapacity_T3mT1','day'});
writetable(T,fullfile(dataDir,'fig5c_scatter_plotdata.csv'));
fprintf('\nDONE fig5.\n');

%% ==================== helpers ====================
function tf=hasTierCap(capDir,DB,lickCol,sess)
  tf=isfile(fullfile(capDir,sprintf('%s_%s_%s_T1_capacity_timecourse.mat',sess,DB,lickCol))) && ...
     isfile(fullfile(capDir,sprintf('%s_%s_%s_T3_capacity_timecourse.mat',sess,DB,lickCol)));
end
function [capT1,capT3,tvec]=loadTierCap(capDir,DB,lickCol,sess)
  S1=load(fullfile(capDir,sprintf('%s_%s_%s_T1_capacity_timecourse.mat',sess,DB,lickCol)),'capacity_rolling','tvec');
  S3=load(fullfile(capDir,sprintf('%s_%s_%s_T3_capacity_timecourse.mat',sess,DB,lickCol)),'capacity_rolling','tvec');
  tvec=S1.tvec(:); capT1=S1.capacity_rolling(:); capT3=S3.capacity_rolling(:);
end
function a=aucVec(aucMap,sess,trials)
  a=nan(numel(trials),1);
  for i=1:numel(trials), k=[sess '|' num2str(trials(i))]; if isKey(aucMap,k), a(i)=aucMap(k); end, end
end
function [tGrid,X,Vstim,auc,keep]=lickMatrix(lickDir,aucMap,DB,sess,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec)
  tGrid=[];X=[];Vstim=[];auc=[];keep=[];
  f=fullfile(lickDir,sprintf('%s_%s_lickTraces.mat',sess,DB)); if ~isfile(f), return; end
  S=load(f); Trial=S.Trial; Vstim=double(S.Visual_Stim); dtv=double(S.Lick_dt); tr=S.lickTrace; n=numel(tr);
  auc=aucVec(aucMap,sess,Trial);
  dtSess=median(dtv(isfinite(dtv)&dtv>0)); if ~isfinite(dtSess)||dtSess<=0, tGrid=[]; return; end
  tGrid=(tPlotRel(1):dtSess:tPlotRel(2))'; nT=numel(tGrid); X=nan(n,nT); keep=false(n,1);
  for i=1:n
    v=tr{i}; dt=dtv(i); if isempty(v)||~isfinite(dt)||dt<=0, continue; end; v=double(v(:));
    if doSmooth&&smoothWinSec>0, v=movmean(v,max(1,round(smoothWinSec/dt)),'Endpoints','shrink'); end
    tRel=(0:numel(v)-1)'*dt-stimOnsetSec; X(i,:)=interp1(tRel,v,tGrid,'linear',NaN); keep(i)=any(isfinite(X(i,:)));
  end
end
function [tGrid,mGoT1,mGoT3,mNoT1,mNoT3]=meanLickTierGoNo(lickDir,aucMap,DB,sess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec)
  mGoT1=[];mGoT3=[];mNoT1=[];mNoT3=[];
  [tGrid,X,Vstim,auc,keep]=lickMatrix(lickDir,aucMap,DB,sess,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec);
  if isempty(tGrid), return; end
  ok=keep&isfinite(auc); if nnz(ok)<10, tGrid=[]; return; end
  isGo=(Vstim==goVal); isNo=~isGo; N=numel(auc);
  [gT1,gT3]=tierMask(auc,ok&isGo); [nT1,nT3]=tierMask(auc,ok&isNo);
  mGoT1=mean(X(gT1,:),1,'omitnan')'; mGoT3=mean(X(gT3,:),1,'omitnan')';
  mNoT1=mean(X(nT1,:),1,'omitnan')'; mNoT3=mean(X(nT3,:),1,'omitnan')';
end
function [tGrid,mAllT1,mAllT3]=meanLickTierAll(lickDir,aucMap,DB,sess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec)
  mAllT1=[];mAllT3=[];
  [tGrid,X,Vstim,auc,keep]=lickMatrix(lickDir,aucMap,DB,sess,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec);
  if isempty(tGrid), return; end
  ok=keep&isfinite(auc); if nnz(ok)<10, tGrid=[]; return; end
  isGo=(Vstim==goVal);
  [gT1,gT3]=tierMask(auc,ok&isGo); [nT1,nT3]=tierMask(auc,ok&~isGo);
  mAllT1=mean(X(gT1|nT1,:),1,'omitnan')'; mAllT3=mean(X(gT3|nT3,:),1,'omitnan')';
end
function [t1,t3]=tierMask(auc,sel)
  t1=false(numel(auc),1); t3=false(numel(auc),1);
  if any(sel), q33=percentileCont(auc(sel),1/3); q67=percentileCont(auc(sel),2/3);
    t1(sel)=auc(sel)<=q33; t3(sel)=auc(sel)>=q67; end
end
function [dCap,dLick,dDay]=collectScatter(capDir,lickDir,aucMap,DB,lickCol,daySess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec,win)
  dCap=[];dLick=[];dDay=[];
  for d=1:6, for si=1:numel(daySess{d})
    sess=char(daySess{d}{si}); if ~hasTierCap(capDir,DB,lickCol,sess), continue; end
    [c1,c3,tvec]=loadTierCap(capDir,DB,lickCol,sess); iC=tvec>=win(1)&tvec<=win(2); if ~any(iC), continue; end
    dc=mean(c3(iC)-c1(iC),'omitnan');
    [tg,l1,l3]=meanLickTierAll(lickDir,aucMap,DB,sess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec);
    if isempty(tg), continue; end; iL=tg>=win(1)&tg<=win(2); if ~any(iL), continue; end
    dl=mean(l3(iL)-l1(iL),'omitnan');
    dCap(end+1,1)=dc; dLick(end+1,1)=dl; dDay(end+1,1)=d; %#ok<AGROW>
  end, end
end
function [c1m,c3m,l1m,l3m,dDay]=collectWinMeans(capDir,lickDir,aucMap,DB,lickCol,daySess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec,win)
  c1m=[];c3m=[];l1m=[];l3m=[];dDay=[];
  for d=1:6, for si=1:numel(daySess{d})
    sess=char(daySess{d}{si}); if ~hasTierCap(capDir,DB,lickCol,sess), continue; end
    [c1,c3,tvec]=loadTierCap(capDir,DB,lickCol,sess); iC=tvec>=win(1)&tvec<=win(2); if ~any(iC), continue; end
    [tg,l1,l3]=meanLickTierAll(lickDir,aucMap,DB,sess,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec);
    if isempty(tg), continue; end; iL=tg>=win(1)&tg<=win(2); if ~any(iL), continue; end
    c1m(end+1,1)=mean(c1(iC),'omitnan'); c3m(end+1,1)=mean(c3(iC),'omitnan'); %#ok<AGROW>
    l1m(end+1,1)=mean(l1(iL),'omitnan'); l3m(end+1,1)=mean(l3(iL),'omitnan'); dDay(end+1,1)=d; %#ok<AGROW>
  end, end
end
function addWindowShade(ax,win,a,rgb)
  yl=get(ax,'YLim'); h=patch(ax,[win(1) win(2) win(2) win(1)],[yl(1) yl(1) yl(2) yl(2)],rgb,'FaceAlpha',a,'EdgeColor','none'); uistack(h,'bottom');
end
function addDiffPatch(ax,t,y1,y2,win,a,rgb)
  t=t(:);y1=y1(:);y2=y2(:); ok=isfinite(t)&isfinite(y1)&isfinite(y2); if nnz(ok)<3, return; end
  tt=t(ok);aa=y1(ok);bb=y2(ok); w=tt>=win(1)&tt<=win(2); if nnz(w)<3, return; end
  tt=tt(w);aa=aa(w);bb=bb(w); h=patch(ax,[tt;flipud(tt)],[aa;flipud(bb)],rgb,'FaceAlpha',a,'EdgeColor','none'); uistack(h,'bottom');
end
function [tRef,y1A,y2A]=alignToReference(tRef,t,y1,y2)
  tRef=tRef(:);t=t(:);y1=y1(:);y2=y2(:);
  if isempty(t)||isempty(y1)||isempty(y2), y1A=nan(size(tRef));y2A=nan(size(tRef)); return; end
  if numel(t)==numel(tRef)&&all(abs(t-tRef)<1e-9), y1A=y1;y2A=y2; return; end
  y1A=interp1(t,y1,tRef,'linear',NaN); y2A=interp1(t,y2,tRef,'linear',NaN);
end
function q=percentileCont(x,p)
  x=x(:); x=x(isfinite(x)); n=numel(x);
  if n==0, q=NaN; return; elseif n==1, q=x(1); return; end
  x=sort(x); if p<=0,q=x(1);return; end; if p>=1,q=x(end);return; end
  r=1+p*(n-1); k=floor(r); dd=r-k; if k>=n, q=x(n); else, q=x(k)+dd*(x(k+1)-x(k)); end
end
function pctLabel(ax,t,y1,y3,win,col)
  % capacity T3-T1: raw difference with symmetric % in parentheses, in the diff-patch colour
  [raw,p]=rawpct(t,y1,y3,win); if isnan(raw), return; end
  yl=ylim(ax);
  text(ax,mean(win),yl(2)-0.06*diff(yl),sprintf('\\Delta_{1-2s} %+.3g (%+.0f%%)',raw,p), ...
    'Color',col,'HorizontalAlignment','center','FontWeight','bold','FontSize',9);
end
function lickLabel(ax,t,gT1,gT3,nT1,nT3,win,col)
  % lick T3-T1, Go and NoGo separately: raw (percentage)
  [rg,pg]=rawpct(t,gT1,gT3,win); [rn,pn]=rawpct(t,nT1,nT3,win);
  yl=ylim(ax);
  text(ax,mean(win),yl(2)-0.02*diff(yl), ...
    {sprintf('Go %+.2f (%+.0f%%)',rg,pg), sprintf('NoGo %+.2f (%+.0f%%)',rn,pn)}, ...
    'Color',col,'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold','FontSize',8);
end
function [raw,p]=rawpct(t,y1,y3,win)
  t=t(:)'; y1=y1(:)'; y3=y3(:)'; w=t>=win(1)&t<=win(2);
  raw=NaN; p=NaN; if ~any(w), return; end
  m1=mean(y1(w),'omitnan'); m3=mean(y3(w),'omitnan'); den=0.5*(abs(m3)+abs(m1));
  raw=m3-m1; if isfinite(den)&&den~=0, p=100*raw/den; end
end
function saveFig(fig,base)
  exportgraphics(fig,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
  exportgraphics(fig,[base '.png'],'Resolution',300,'BackgroundColor','white');
  savefig(fig,[base '.fig']);
  fprintf('saved %s\n',base);
end
