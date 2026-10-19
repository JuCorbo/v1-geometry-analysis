%% make_fig5_pctBoxes_combined.m  —  Figure 5 main panel d (ENHANCED)
% One summary boxplot: each behavioral variable (Licking, Locomotion, Pupil pre-stim)
% differs strongly T3-T1, while discrimination Capacity barely moves for any split;
% Capacity Trained-Naive (random thirds) is the training-effect reference.
clear; clc; close all;
codeDir=fileparts(mfilename('fullpath')); anaDir=fileparts(codeDir); figDir=fileparts(fileparts(anaDir));
capDir=fullfile(figDir,'Figure4','manifoldCapacity','cache','capacityTimecourse'); lickDir=fullfile(anaDir,'cache','lickTracesFull');
aucCsv=fullfile(anaDir,'cache','LickAUCz_perTrial.csv');
locoCache=fullfile(anaDir,'cache','loco'); pupCache=fullfile(anaDir,'cache','pupil');
sharedTr=fullfile(figDir,'_sharedWindowedData','data','trained'); sharedNv=fullfile(figDir,'_sharedWindowedData','data','naive');
locoQCdir=fullfile(anaDir,'cache','locoQC');
pupQCdir=fullfile(anaDir,'cache','pupilQC');
panelDir=fullfile(anaDir,'panels','main'); statsDir=fullfile(anaDir,'stats');
DB='BatuhanAngDist'; dbNaive='BatuhanADNaive'; lickCol='LickAUCz_1_2'; goVal=45;
days={'D1','D2','D3','D4','D5','D6'}; quantWin=[1 2]; stimOnsetSec=1.5; tPlotRel=[-1 4]; doSmooth=true; smoothWinSec=0.15;
dt=0.0646; fs=1/dt;

daySess=cell(6,1); naiveDaySess=cell(6,1);
for d=1:6
  ff=dir(fullfile(sharedTr,sprintf('D%d',d),'*_win.mat')); daySess{d}=erase(erase({ff.name},[DB '_']),'_win.mat');
  gg=dir(fullfile(sharedNv,sprintf('D%d',d),'*_win.mat')); naiveDaySess{d}=erase(erase({gg.name},[dbNaive '_']),'_win.mat');
end
Q=readtable(aucCsv); aucMap=containers.Map('KeyType','char','ValueType','double');
for i=1:height(Q), aucMap([char(string(Q.Experiment(i))) '|' num2str(Q.Trial(i))])=Q.LickAUCz_1_2(i); end

%% ---- LICK tertiles: %Delta licking + capacity (all days) ----
pctLick=[]; pctCapL=[]; LkW1=[];LkW3=[];CpL1=[];CpL3=[];
for d=1:6, for si=1:numel(daySess{d})
  se=char(daySess{d}{si}); if ~hasTier(capDir,DB,lickCol,se), continue; end
  [c1,c3,tv]=loadTier(capDir,DB,lickCol,se); iC=tv>=quantWin(1)&tv<=quantWin(2);
  [tg,l1,l3]=meanLickTierAll(lickDir,aucMap,DB,se,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec);
  if isempty(tg), continue; end; iL=tg>=quantWin(1)&tg<=quantWin(2);
  cap1=mean(c1(iC),'omitnan'); cap3=mean(c3(iC),'omitnan'); lk1=mean(l1(iL),'omitnan'); lk3=mean(l3(iL),'omitnan');
  pctLick(end+1,1)=100*(lk3-lk1)/(0.5*(abs(lk3)+abs(lk1))); %#ok<SAGROW>
  pctCapL(end+1,1)=100*(cap3-cap1)/(0.5*(abs(cap3)+abs(cap1))); %#ok<SAGROW>
  LkW1(end+1,1)=lk1;LkW3(end+1,1)=lk3;CpL1(end+1,1)=cap1;CpL3(end+1,1)=cap3; %#ok<SAGROW>
end, end

%% ---- Trained-Naive (random thirds) capacity %Delta ----
pctTrNv=[];
for d=1:6
  nm=[]; for si=1:numel(naiveDaySess{d}), s=char(naiveDaySess{d}{si}); f=fullfile(capDir,sprintf('%s_%s_randThirdAvg_capacity_timecourse.mat',s,dbNaive)); if ~isfile(f), continue; end; S=load(f); iw=S.tvec(:)>=quantWin(1)&S.tvec(:)<=quantWin(2); nm(end+1,1)=mean(S.capacity_rolling_randThirdAvg(iw),'omitnan'); end
  ndm=mean(nm,'omitnan');
  for si=1:numel(daySess{d}), s=char(daySess{d}{si}); f=fullfile(capDir,sprintf('%s_%s_randThirdAvg_capacity_timecourse.mat',s,DB)); if ~isfile(f), continue; end; S=load(f); iw=S.tvec(:)>=quantWin(1)&S.tvec(:)<=quantWin(2); tw=mean(S.capacity_rolling_randThirdAvg(iw),'omitnan'); den=0.5*(abs(tw)+abs(ndm)); if den==0, continue; end; pctTrNv(end+1,1)=100*(tw-ndm)/den; end
end

%% ---- LOCO tertiles: %Delta loco + capacity ----
QL=readtable(fullfile(locoQCdir,sprintf('locoQC_%s.csv',DB))); wl=QL(QL.pass_post==1 & QL.sessStd>=3,{'Experiment','BehavCond'});
excl="F19_42121_"; StL=readtable(fullfile(locoQCdir,sprintf('locoStereotypy_%s.csv',DB))); excl=union(excl,string(StL.sess(StL.R2mean>0.30))); wl=wl(~ismember(string(wl.Experiment),excl),:);
bLoco=[]; bCapLoco=[]; LoW1=[];LoW3=[];CpO1=[];CpO3=[]; winSec=[1 2]; minFrames=60; Ncommon=90;
for r=1:height(wl)
  se=char(wl.Experiment(r)); day=char(wl.BehavCond(r));
  c1=fullfile(locoCache,sprintf('%s_%s_loco_post_T1_capacity_timecourse.mat',se,DB)); c3=fullfile(locoCache,sprintf('%s_%s_loco_post_T3_capacity_timecourse.mat',se,DB)); prep=fullfile(locoCache,sprintf('%s_%s_prep.mat',se,DB));
  if ~all(isfile({c1,c3,prep})), continue; end
  C1=load(c1); C3=load(c3); tv=C1.tvec(:)'; cap1=C1.capacity_rolling(:)'; cap3=C3.capacity_rolling(:)'; P=load(prep,'T1','T3','isGo');
  [Mn,xval,isGo,valid]=fetchLocoNAS(locoCache,se,DB,goVal,dt,fs,stimOnsetSec,1,minFrames,P.isGo,Ncommon);
  if isempty(Mn), continue; end
  wlm=xval>=winSec(1)&xval<=winSec(2); perTrial=mean(Mn(:,wlm),2); l1=mean(perTrial(P.T1(valid(P.T1))),'omitnan'); l3=mean(perTrial(P.T3(valid(P.T3))),'omitnan');
  bLoco(end+1,1)=100*(l3-l1)/(0.5*(abs(l3)+abs(l1))); wi=tv>=winSec(1)&tv<=winSec(2); den=0.5*(cap3+cap1); pv=100*(cap3-cap1)./den; pv(~isfinite(pv)|den==0)=NaN; bCapLoco(end+1,1)=mean(pv(wi),'omitnan');
  LoW1(end+1,1)=l1;LoW3(end+1,1)=l3;CpO1(end+1,1)=mean(cap1(wi),'omitnan');CpO3(end+1,1)=mean(cap3(wi),'omitnan');
end

%% ---- PUPIL pre-stim (baseline split): %Delta pupil + capacity ----
QP=readtable(fullfile(pupQCdir,sprintf('pupilQC_%s.csv',DB))); wp=QP(QP.pass_baseline==1,{'Experiment','BehavCond'}); wp=wp(~ismember(string(wp.Experiment),"F10_42521_"),:);
qc=struct('maxNaNFracWindow',0.10,'jumpMADk',8,'maxJumpFracTrial',0.02,'minStdRaw',1e-6); winSecP=[-1.5 0];
bPup=[]; bCapPup=[]; PuW1=[];PuW3=[];CpP1=[];CpP3=[];
for r=1:height(wp)
  se=char(wp.Experiment(r));
  c1=fullfile(pupCache,sprintf('%s_%s_pupil_baseline_T1_capacity_timecourse.mat',se,DB)); c3=fullfile(pupCache,sprintf('%s_%s_pupil_baseline_T3_capacity_timecourse.mat',se,DB));
  if ~(isfile(c1)&&isfile(c3)), continue; end
  [Mn,xval,isGo,valid,metric,usable]=computePupilNAS(pupCache,se,DB,goVal,dt,fs,stimOnsetSec,2,minFrames,Ncommon,winSecP,qc);
  if isempty(Mn), continue; end
  zu=metric(usable); z=(metric-mean(zu,'omitnan'))/std(zu,'omitnan'); [T1,T3]=withinClassTertiles(z,usable,isGo);
  C1=load(c1); C3=load(c3); tv=C1.tvec(:)'; cap1=C1.capacity_rolling(:)'; cap3=C3.capacity_rolling(:)';
  wlm=xval>=winSecP(1)&xval<=winSecP(2); perTrial=mean(Mn(:,wlm),2,'omitnan'); p1=mean(perTrial(T1(valid(T1))),'omitnan'); p3=mean(perTrial(T3(valid(T3))),'omitnan');
  bPup(end+1,1)=100*(p3-p1)/(0.5*(abs(p3)+abs(p1))); wi=tv>=quantWin(1)&tv<=quantWin(2); den=0.5*(cap3+cap1); pv=100*(cap3-cap1)./den; pv(~isfinite(pv)|den==0)=NaN; bCapPup(end+1,1)=mean(pv(wi),'omitnan');
  PuW1(end+1,1)=p1;PuW3(end+1,1)=p3;CpP1(end+1,1)=mean(cap1(wi),'omitnan');CpP3(end+1,1)=mean(cap3(wi),'omitnan');
end

%% ---- combined boxplot ----
labels={'Lick','Locomotion','Pupil(pre)','Cap|lick','Cap|loco','Cap|pupil','Cap Tr-Nv'};
data={pctLick,bLoco,bPup,pctCapL,bCapLoco,bCapPup,pctTrNv};
colBeh=[0.15 0.6 0.2]; colCap=[0.25 0.25 0.25]; colTr=[0.2 0.4 0.8];
cols={colBeh,colBeh,colBeh,colCap,colCap,colCap,colTr};
fig=figure('Color','w','Position',[60 60 1150 620]); ax=axes(fig); hold(ax,'on');
for i=1:numel(data)
  v=data{i}; v=v(isfinite(v));
  boxchart(ax,repmat(i,numel(v),1),v,'BoxWidth',0.55,'BoxFaceColor',cols{i},'MarkerStyle','none','WhiskerLineColor',cols{i});
  scatter(ax,i+0.14*(rand(numel(v),1)-0.5),v,22,cols{i},'filled','MarkerFaceAlpha',0.5);
  text(ax,i,max(v)+6,sprintf('%.0f%%\nn=%d',median(v),numel(v)),'HorizontalAlignment','center','FontSize',9);
end
yline(ax,0,'k-'); xline(ax,3.5,'k:','HandleVisibility','off'); xline(ax,6.5,'k:','HandleVisibility','off');
set(ax,'XTick',1:numel(labels),'XTickLabel',labels,'FontSize',11); xlim(ax,[0.4 numel(labels)+0.6]);
ylabel(ax,'T3-T1 symmetric % change (behavior in its split window)');
title(ax,'Behavioral variables differ strongly between tertiles; discrimination capacity does not (training does)','FontWeight','bold');
text(ax,2,-1,'behavioral driver','HorizontalAlignment','center','Color',colBeh,'FontWeight','bold','VerticalAlignment','top');
text(ax,5,-1,'capacity (per split)','HorizontalAlignment','center','Color',colCap,'FontWeight','bold','VerticalAlignment','top');
box(ax,'off'); grid(ax,'on');
if ~exist(panelDir,'dir'), mkdir(panelDir); end
exportgraphics(fig,fullfile(panelDir,'fig5d_pctBoxes_combined.png'),'Resolution',220,'BackgroundColor','white');
exportgraphics(fig,fullfile(panelDir,'fig5d_pctBoxes_combined.pdf'),'ContentType','vector','BackgroundColor','white');
if ~exist(statsDir,'dir'), mkdir(statsDir); end
fid=fopen(fullfile(statsDir,'fig5d_pctBoxes_combined_stats.txt'),'w');
for i=1:numel(labels), v=data{i}; v=v(isfinite(v)); fprintf(fid,'%-12s median=%.1f%% IQR[%.1f, %.1f] n=%d signrank-p=%.3g\n',labels{i},median(v),prctile(v,25),prctile(v,75),numel(v),signrank(v)); end
fprintf(fid,'\nHedges g_av (paired T3 vs T1, raw [1-2]s window means):\n');
fprintf(fid,'  Lick  g=%.2f | Cap|lick  g=%.2f\n', gav(LkW1,LkW3), gav(CpL1,CpL3));
fprintf(fid,'  Loco  g=%.2f | Cap|loco  g=%.2f\n', gav(LoW1,LoW3), gav(CpO1,CpO3));
fprintf(fid,'  Pupil g=%.2f | Cap|pupil g=%.2f\n', gav(PuW1,PuW3), gav(CpP1,CpP3));
fclose(fid); type(fullfile(statsDir,'fig5d_pctBoxes_combined_stats.txt'));
fprintf('saved combined pct boxes\n');

%% ---- helpers (reused) ----
function g=gav(v1,v3)
  ok=isfinite(v1)&isfinite(v3); v1=v1(ok); v3=v3(ok); N=numel(v1);
  if N<3, g=NaN; return; end
  sdav=sqrt((std(v1)^2+std(v3)^2)/2); if sdav==0, g=NaN; return; end
  g=(1-3/(4*N-9))*(mean(v3-v1)/sdav);
end
function tf=hasTier(capDir,DB,lc,se), tf=isfile(fullfile(capDir,sprintf('%s_%s_%s_T1_capacity_timecourse.mat',se,DB,lc)))&&isfile(fullfile(capDir,sprintf('%s_%s_%s_T3_capacity_timecourse.mat',se,DB,lc))); end
function [c1,c3,tv]=loadTier(capDir,DB,lc,se), S1=load(fullfile(capDir,sprintf('%s_%s_%s_T1_capacity_timecourse.mat',se,DB,lc)),'capacity_rolling','tvec'); S3=load(fullfile(capDir,sprintf('%s_%s_%s_T3_capacity_timecourse.mat',se,DB,lc)),'capacity_rolling','tvec'); tv=S1.tvec(:)'; c1=S1.capacity_rolling(:)'; c3=S3.capacity_rolling(:)'; end
function a=aucVec(aucMap,se,tr), a=nan(numel(tr),1); for i=1:numel(tr), k=[se '|' num2str(tr(i))]; if isKey(aucMap,k), a(i)=aucMap(k); end, end, end
function [tg,mAllT1,mAllT3]=meanLickTierAll(lickDir,aucMap,DB,se,goVal,stimOnsetSec,tPlotRel,doSmooth,smoothWinSec)
  tg=[];mAllT1=[];mAllT3=[]; f=fullfile(lickDir,sprintf('%s_%s_lickTraces.mat',se,DB)); if ~isfile(f), return; end
  S=load(f); Trial=S.Trial; Vs=double(S.Visual_Stim); dtv=double(S.Lick_dt); tr=S.lickTrace; n=numel(tr); auc=aucVec(aucMap,se,Trial);
  dtS=median(dtv(isfinite(dtv)&dtv>0)); if ~(dtS>0), tg=[]; return; end
  tg=(tPlotRel(1):dtS:tPlotRel(2))'; nT=numel(tg); X=nan(n,nT); keep=false(n,1);
  for i=1:n, v=tr{i}; if isempty(v), continue; end, v=double(v(:)); if doSmooth, v=movmean(v,max(1,round(smoothWinSec/dtS)),'Endpoints','shrink'); end, tRel=(0:numel(v)-1)'*dtS-stimOnsetSec; X(i,:)=interp1(tRel,v,tg,'linear',NaN); keep(i)=any(isfinite(X(i,:))); end
  ok=keep&isfinite(auc); if nnz(ok)<10, tg=[]; return; end
  isGo=Vs==goVal; [g1,g3]=tierMask(auc,ok&isGo); [n1,n3]=tierMask(auc,ok&~isGo);
  mAllT1=mean(X(g1|n1,:),1,'omitnan')'; mAllT3=mean(X(g3|n3,:),1,'omitnan')';
end
function [t1,t3]=tierMask(auc,sel), t1=false(numel(auc),1); t3=false(numel(auc),1); if any(sel), q33=pctl(auc(sel),1/3); q67=pctl(auc(sel),2/3); t1(sel)=auc(sel)<=q33; t3(sel)=auc(sel)>=q67; end, end
function q=pctl(x,p), x=x(:); x=x(isfinite(x)); n=numel(x); if n==0,q=NaN;return;elseif n==1,q=x(1);return;end, x=sort(x); if p<=0,q=x(1);return;end; if p>=1,q=x(end);return;end, r=1+p*(n-1); k=floor(r); d=r-k; if k>=n,q=x(n);else,q=x(k)+d*(x(k+1)-x(k));end, end
function [T1,T3]=withinClassTertiles(z,usable,isGo), T1=[];T3=[]; for cls=[true false], idx=find(usable&(isGo==cls)); if numel(idx)<3, continue; end, vals=z(idx); q33=quantile(vals,1/3); q67=quantile(vals,2/3); T1=[T1; idx(vals<=q33)]; T3=[T3; idx(vals>=q67)]; end, T1=sort(T1); T3=sort(T3); end
function [Mn,xval,isGo,valid]=fetchLocoNAS(locoDir,se,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,prepIsGo,Ncommon)
  Mn=[];xval=[];isGo=[];valid=[]; f=fullfile(locoDir,sprintf('%s_%s_locoTraces.mat',se,DB)); if ~isfile(f), return; end
  S=load(f); uTr=double(S.Trial); isGo=double(S.Visual_Stim)==goVal;
  if ~isempty(prepIsGo)&&(numel(prepIsGo)~=numel(uTr)||~isequal(logical(prepIsGo(:)),isGo(:))), return; end
  rb=S.locoTrace; lens=cellfun(@numel,rb); val=lens>=minFrames & ~cellfun(@(v) isempty(v)||all(~isfinite(v)),rb); if nnz(val)<6, return; end
  nF=min(lens(val)); gi=find(val(:))'; M=NaN(numel(uTr),nF); for a=gi, v=double(rb{a}(:))'; v=v(1:nF); v(~isfinite(v))=0; M(a,:)=v; end
  Mlp=NaN(size(M)); Mlp(gi,:)=lowpass(M(gi,:)',lowpassHz,fs)'; sd=std(Mlp(gi,:),0,'all','omitnan'); if sd>0, Mn=Mlp/sd; else, Mn=Mlp; end
  if size(Mn,2)>=Ncommon, Mn=Mn(:,1:Ncommon); else, Mn(:,end+1:Ncommon)=NaN; end
  xval=(0:Ncommon-1)*dt-stimOnsetSec; valid=false(numel(uTr),1); valid(gi)=true;
end
function [Mn,xval,isGo,valid,metric,usable]=computePupilNAS(pupDir,se,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,Ncommon,winSec,qc)
  Mn=[];xval=[];isGo=[];valid=[];metric=[];usable=[]; f=fullfile(pupDir,sprintf('%s_%s_pupilTraces.mat',se,DB)); if ~isfile(f), return; end
  S=load(f); uTr=double(S.Trial); isGo=double(S.Visual_Stim)==goVal; rb=S.pupilTrace; n=numel(uTr);
  LP=cell(n,1); metric=NaN(n,1); usable=false(n,1); lens=zeros(n,1);
  for a=1:n
    v=rb{a}; if isempty(v)||all(isnan(v)), continue; end, v=double(v(:)'); L=numel(v); lens(a)=L; if L<minFrames, continue; end
    tRel=(0:L-1)*dt-stimOnsetSec; idx=tRel>=winSec(1)&tRel<=winSec(2); vi=v; if any(isnan(vi)), g=find(~isnan(vi)); if numel(g)<3, continue; end, vi=interp1(g,vi(g),1:L,'linear','extrap'); end
    Wn=lowpassHz/(fs/2); vLP=vi; if Wn<1, [bb,aa]=butter(2,Wn,'low'); vLP=filtfilt(bb,aa,vi); end, LP{a}=vLP;
    if any(idx), cov=mean(~isnan(v(idx))); dd=diff(vi); okjump=true; if numel(dd)>=5, m=mad(dd,1); jf=(m>0)*mean(abs(dd-median(dd))>qc.jumpMADk*m); okjump=jf<=qc.maxJumpFracTrial; end
      if cov>=(1-qc.maxNaNFracWindow)&&std(vi(idx),0,'omitnan')>qc.minStdRaw&&okjump, metric(a)=mean(vLP(idx),'omitnan'); usable(a)=true; end, end
  end
  gi=find(cellfun(@(x)~isempty(x),LP))'; if numel(gi)<6, return; end
  nF=min(lens(gi)); M=NaN(n,nF); for a=gi, M(a,:)=LP{a}(1:nF); end, sd=std(M(gi,:),0,'all','omitnan'); if sd>0, M=M/sd; end
  if size(M,2)>=Ncommon, Mn=M(:,1:Ncommon); else, Mn=[M NaN(n,Ncommon-size(M,2))]; end
  xval=(0:Ncommon-1)*dt-stimOnsetSec; valid=false(n,1); valid(gi)=true;
end
