%% make_fig5_supp_loco.m  —  Figure S10: locomotion splits on capacity
% Data:
%   - loco T1/T3/EQ capacity + locoPrep splits : cache/loco/
%   - all-trials capacity                       : Figure4 capacityTimecourse cache
%   - locomotion traces                         : cache/loco/<sess>_locoTraces.mat
%   - working set (QC)                          : cache/locoQC
clear; clc; close all;
codeDir=fileparts(mfilename('fullpath')); anaDir=fileparts(codeDir);
locoCache=fullfile(anaDir,'cache','loco');
allCapDir=fullfile(fileparts(fileparts(anaDir)),'Figure4','manifoldCapacity','cache','capacityTimecourse');
qcDir=fullfile(anaDir,'cache','locoQC');
panelDir=fullfile(anaDir,'panels','supp'); statsDir=fullfile(anaDir,'stats');

DB='BatuhanAngDist'; goVal=45; days={'D1','D2','D3','D4','D5','D6'}; noGoVal=[135,90,75,70,65,60];
dt=0.0646; fs=1/dt; stimOnsetSec=1.5; lowpassHz=1; minFrames=60;
winSec=[1 2]; xlimSec=[-1 4]; Ncommon=90; exampleSess='C05_121020'; exampleDay='D1';
colT1=[0.85 0.15 0.15]; colT3=[0.10 0.30 0.85]; colEQ=[0.50 0.15 0.60]; colAll=[0.15 0.15 0.15];
colGo=[0.20 0.45 0.80]; colNo=[0.85 0.30 0.25]; lsGo='--'; lsNo=':';

QCt=readtable(fullfile(qcDir,sprintf('locoQC_%s.csv',DB)));
work=QCt(QCt.pass_post==1 & QCt.sessStd>=3, {'Experiment','BehavCond'});
excl="F19_42121_"; St=readtable(fullfile(qcDir,sprintf('locoStereotypy_%s.csv',DB)));
excl=union(excl,string(St.sess(St.R2mean>0.30))); work=work(~ismember(string(work.Experiment),excl),:);
fprintf('Working set: %d\n',height(work));

capT1=cell(6,1);capT3=cell(6,1);LT1Go=cell(6,1);LT1No=cell(6,1);LT3Go=cell(6,1);LT3No=cell(6,1);
capEQ=cell(6,1);capAll=cell(6,1);EQGo=cell(6,1);EQNo=cell(6,1); tvec=[]; bLoco=[];bCap=[];bAllEq=[]; skipped={};

% example (panel A)
[exGo,exNo,exGoM,exNoM,exX]=exLoco(locoCache,exampleSess,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,Ncommon);

for r=1:height(work)
  se=char(work.Experiment(r)); day=char(work.BehavCond(r)); d=find(strcmp(days,day));
  c1=fullfile(locoCache,sprintf('%s_%s_loco_post_T1_capacity_timecourse.mat',se,DB));
  c3=fullfile(locoCache,sprintf('%s_%s_loco_post_T3_capacity_timecourse.mat',se,DB));
  ceq=fullfile(locoCache,sprintf('%s_%s_loco_post_EQ_capacity_timecourse.mat',se,DB));
  cal=fullfile(allCapDir,sprintf('%s_%s_capacity_timecourse.mat',se,DB));
  prep=fullfile(locoCache,sprintf('%s_%s_prep.mat',se,DB));
  if ~all(isfile({c1,c3,ceq,cal,prep})), skipped{end+1}=[se ' (missing cache)']; continue; end
  C1=load(c1);C3=load(c3);CE=load(ceq);CA=load(cal); tvec=C1.tvec(:)';
  wi=tvec>=winSec(1)&tvec<=winSec(2); bi=tvec<0;
  cap1=C1.capacity_rolling(:)'; cap3=C3.capacity_rolling(:)';
  P=load(prep,'T1','T3','EQ','isGo');
  [Mn,xval,isGo,valid]=fetchLocoNAS(locoCache,se,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,P.isGo,Ncommon);
  if isempty(Mn), skipped{end+1}=[se ' (no loco)']; continue; end
  capT1{d}=[capT1{d};cap1]; capT3{d}=[capT3{d};cap3];
  capEQ{d}=[capEQ{d}; CE.capacity_rolling(:)'-mean(CE.capacity_rolling(bi),'omitnan')];
  capAll{d}=[capAll{d}; CA.capacity_rolling(:)'-mean(CA.capacity_rolling(bi),'omitnan')];
  sub=@(idx) mean(Mn(idx(valid(idx)),:),1,'omitnan');
  LT1Go{d}=[LT1Go{d};sub(P.T1(isGo(P.T1)))]; LT1No{d}=[LT1No{d};sub(P.T1(~isGo(P.T1)))];
  LT3Go{d}=[LT3Go{d};sub(P.T3(isGo(P.T3)))]; LT3No{d}=[LT3No{d};sub(P.T3(~isGo(P.T3)))];
  EQGo{d}=[EQGo{d};sub(P.EQ(isGo(P.EQ)))]; EQNo{d}=[EQNo{d};sub(P.EQ(~isGo(P.EQ)))];
  wl=xval>=winSec(1)&xval<=winSec(2); perTrial=mean(Mn(:,wl),2);
  l1=mean(perTrial(P.T1(valid(P.T1))),'omitnan'); l3=mean(perTrial(P.T3(valid(P.T3))),'omitnan');
  bLoco(end+1,1)=100*(l3-l1)/(0.5*(l3+l1));
  den=0.5*(cap3+cap1); pv=100*(cap3-cap1)./den; pv(~isfinite(pv)|den==0)=NaN; bCap(end+1,1)=mean(pv(wi),'omitnan');
  evEQ=mean(capEQ{d}(end,wi),'omitnan'); evAL=mean(capAll{d}(end,wi),'omitnan'); bAllEq(end+1,1)=100*(evAL-evEQ)/(0.5*(evAL+evEQ));
end
fprintf('box n=%d | loco med=%.0f%% cap med=%.1f%% allEq med=%.1f%% | skipped: %s\n',numel(bLoco),median(bLoco),median(bCap),median(bAllEq,'omitnan'),strjoin(skipped,'; '));

%% figure (verbatim layout)
fig=figure('Color','w','Position',[20 20 1040 2300]);
TL=tiledlayout(fig,6,1,'TileSpacing','compact','Padding','compact');
title(TL,sprintf('%s | Locomotion is not a confound for discrimination capacity (n=%d sessions)',DB,numel(bLoco)),'Interpreter','none','FontWeight','bold','FontSize',13);
tTop=tiledlayout(TL,1,2,'TileSpacing','compact','Padding','compact'); tTop.Layout.Tile=1; tTop.Layout.TileSpan=[2 1];
axA=nexttile(tTop); hold(axA,'on');
plotSingles(axA,exX,exGo,colGo); plotSingles(axA,exX,exNo,colNo);
plot(axA,exX,exGoM,'-','Color',colGo,'LineWidth',2.6); plot(axA,exX,exNoM,'-','Color',colNo,'LineWidth',2.6);
xline(axA,0,'k:'); xline(axA,3,'k:'); xlim(axA,xlimSec); xlabel(axA,'Time from stim (s)'); ylabel(axA,'Locomotion (AU, /sessStd)');
title(axA,sprintf('(A) Example %s: single trials + mean (Go blue > NoGo red)',exampleSess),'Interpreter','none'); box(axA,'off'); grid(axA,'on'); axis(axA,'square');
axD=nexttile(tTop); hold(axD,'on');
cats=categorical({'Loco T3-T1','Cap T3-T1','Cap ALL-EQ'},{'Loco T3-T1','Cap T3-T1','Cap ALL-EQ'}); n=numel(bLoco);
X=[repmat(cats(1),n,1);repmat(cats(2),n,1);repmat(cats(3),sum(isfinite(bAllEq)),1)]; Y=[bLoco;bCap;bAllEq(isfinite(bAllEq))];
b=boxchart(axD,X,Y,'BoxWidth',0.5); b.BoxFaceColor=[0.85 0.85 0.85]; b.WhiskerLineColor=[0 0 0]; b.MarkerStyle='none'; rng(1);
sc=@(xc,v) scatter(axD,repmat(double(xc),numel(v),1)+0.12*(rand(numel(v),1)-0.5),v,26,[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.6);
sc(cats(1),bLoco); sc(cats(2),bCap); sc(cats(3),bAllEq(isfinite(bAllEq)));
yline(axD,0,'k-'); ylabel(axD,'T3-T1 symmetric % (1-2 s)');
title(axD,sprintf('(D) Summary: loco med %.0f%% (p=%.1g) | cap %.1f%% (p=%.2g) | ALL-EQ %.1f%%',median(bLoco),signrank(bLoco),median(bCap),signrank(bCap),median(bAllEq,'omitnan')),'Interpreter','tex');
grid(axD,'on'); box(axD,'off'); axis(axD,'square');
wrapB=tiledlayout(TL,1,3,'TileSpacing','none','Padding','none'); wrapB.Layout.Tile=3; wrapB.Layout.TileSpan=[2 1];
tB=tiledlayout(wrapB,2,3,'TileSpacing','compact','Padding','compact'); tB.Layout.Tile=1; tB.Layout.TileSpan=[1 2];   % B panels ~2/3 width (less wide)
title(tB,'(B) T1 (red, low-loco) vs T3 (blue, high-loco): capacity (left, solid) & locomotion (right, Go --/ NoGo :)','FontWeight','bold','FontSize',9);
for d=1:6
  ax=nexttile(tB); hold(ax,'on'); yyaxis(ax,'left');
  for rr=1:size(capT1{d},1), p=plot(ax,tvec,capT1{d}(rr,:),'-','LineWidth',0.5); p.Color=[colT1 0.18]; p.HandleVisibility='off'; end
  for rr=1:size(capT3{d},1), p=plot(ax,tvec,capT3{d}(rr,:),'-','LineWidth',0.5); p.Color=[colT3 0.18]; p.HandleVisibility='off'; end
  plot(ax,tvec,mean(capT1{d},1,'omitnan'),'-','Color',colT1,'LineWidth',2.2); plot(ax,tvec,mean(capT3{d},1,'omitnan'),'-','Color',colT3,'LineWidth',2.2);
  ax.YColor=[0.2 0.2 0.2]; if mod(d,3)==1, ylabel(ax,'Capacity'); end
  yyaxis(ax,'right');
  [g1,xg]=avgAl(LT1Go{d},dt,stimOnsetSec); [n1,~]=avgAl(LT1No{d},dt,stimOnsetSec); [g3,~]=avgAl(LT3Go{d},dt,stimOnsetSec); [n3,~]=avgAl(LT3No{d},dt,stimOnsetSec);
  plot(ax,xg,g1,lsGo,'Color',colT1,'LineWidth',1.4); plot(ax,xg,g3,lsGo,'Color',colT3,'LineWidth',1.4);
  plot(ax,xg,n1,lsNo,'Color',colT1,'LineWidth',1.4); plot(ax,xg,n3,lsNo,'Color',colT3,'LineWidth',1.4);
  ax.YColor=[0.2 0.2 0.2]; if mod(d,3)==0, ylabel(ax,'Loco (AU)'); end
  xline(ax,0,'k:'); xline(ax,3,'k:'); xlim(ax,xlimSec); grid(ax,'on'); box(ax,'off'); title(ax,sprintf('%s | NoGo=%g',days{d},noGoVal(d))); if d==6, xlabel(ax,'Time from stim (s)'); end
end
wrapC=tiledlayout(TL,1,3,'TileSpacing','none','Padding','none'); wrapC.Layout.Tile=5; wrapC.Layout.TileSpan=[2 1];
tC=tiledlayout(wrapC,2,3,'TileSpacing','compact','Padding','compact'); tC.Layout.Tile=1; tC.Layout.TileSpan=[1 2];   % C panels ~2/3 width (less wide)
title(tC,'(C) EQ split: baseline-aligned capacity EQ (purple) vs ALL trials (black, left) + EQ locomotion (Go --/ NoGo :, right)','FontWeight','bold','FontSize',9);
for d=1:6
  ax=nexttile(tC); hold(ax,'on'); yyaxis(ax,'left');
  muE=mean(capEQ{d},1,'omitnan'); seE=std(capEQ{d},0,1,'omitnan')./sqrt(size(capEQ{d},1)); muA=mean(capAll{d},1,'omitnan'); seA=std(capAll{d},0,1,'omitnan')./sqrt(size(capAll{d},1));
  fill(ax,[tvec fliplr(tvec)],[muE+seE fliplr(muE-seE)],colEQ,'EdgeColor','none','FaceAlpha',0.12,'HandleVisibility','off');
  fill(ax,[tvec fliplr(tvec)],[muA+seA fliplr(muA-seA)],colAll,'EdgeColor','none','FaceAlpha',0.12,'HandleVisibility','off');
  plot(ax,tvec,muE,'-','Color',colEQ,'LineWidth',2.2); plot(ax,tvec,muA,'-','Color',colAll,'LineWidth',2.2);
  yline(ax,0,'k-'); ax.YColor=[0.2 0.2 0.2]; if mod(d,3)==1, ylabel(ax,'\Delta Capacity'); end
  yyaxis(ax,'right'); [g,xg]=avgAl(EQGo{d},dt,stimOnsetSec); [nn,~]=avgAl(EQNo{d},dt,stimOnsetSec);
  plot(ax,xg,g,lsGo,'Color',colGo,'LineWidth',1.4); plot(ax,xg,nn,lsNo,'Color',colNo,'LineWidth',1.4);
  ax.YColor=[0.2 0.2 0.2]; if mod(d,3)==0, ylabel(ax,'EQ loco'); end
  xline(ax,0,'k:'); xline(ax,3,'k:'); xlim(ax,xlimSec); grid(ax,'on'); box(ax,'off'); title(ax,sprintf('%s | NoGo=%g',days{d},noGoVal(d))); if d==6, xlabel(ax,'Time from stim (s)'); end
end
if ~exist(panelDir,'dir'), mkdir(panelDir); end
exportgraphics(fig,fullfile(panelDir,sprintf('fig5_supp_loco_%s.png',DB)),'Resolution',180,'BackgroundColor','white');
exportgraphics(fig,fullfile(panelDir,sprintf('fig5_supp_loco_%s.pdf',DB)),'ContentType','vector','BackgroundColor','white');
if ~exist(statsDir,'dir'), mkdir(statsDir); end
fid=fopen(fullfile(statsDir,'fig5_supp_loco_stats.txt'),'w');
fprintf(fid,'Loco supp: n=%d | loco T3-T1 med=%.0f%% (signrank p=%.2g); cap T3-T1 med=%.1f%% (p=%.2g); cap ALL-EQ med=%.1f%%\n',numel(bLoco),median(bLoco),signrank(bLoco),median(bCap),signrank(bCap),median(bAllEq,'omitnan'));
fclose(fid); fprintf('saved loco supp figure + stats\n');

%% helpers
function plotSingles(ax,x,M,c), for k=1:size(M,1), y=M(k,:); if all(isnan(y)), continue; end, plot(ax,x,y,'Color',[c 0.05],'LineWidth',0.3,'HandleVisibility','off'); end, end
function [mu,xg]=avgAl(traces,dt,stimOnsetSec), if isempty(traces), mu=[]; xg=[]; return; end, maxF=size(traces,2); mu=mean(traces,1,'omitnan'); xg=(0:maxF-1)*dt-stimOnsetSec; end
function [exGo,exNo,exGoM,exNoM,xval]=exLoco(locoDir,se,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,Ncommon)
  [Mn,xval,isGo,valid]=fetchLocoNAS(locoDir,se,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,[],Ncommon);
  if isempty(Mn), exGo=[];exNo=[];exGoM=[];exNoM=[];xval=[]; return; end
  exGo=Mn(valid&isGo,:); exNo=Mn(valid&~isGo,:); exGoM=mean(exGo,1,'omitnan'); exNoM=mean(exNo,1,'omitnan');
end
function [Mn,xval,isGo,valid]=fetchLocoNAS(locoDir,se,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,prepIsGo,Ncommon)
  Mn=[];xval=[];isGo=[];valid=[];
  f=fullfile(locoDir,sprintf('%s_%s_locoTraces.mat',se,DB)); if ~isfile(f), return; end
  S=load(f); uTr=double(S.Trial); isGo=double(S.Visual_Stim)==goVal;   % prep-sorted by Trial
  if ~isempty(prepIsGo) && (numel(prepIsGo)~=numel(uTr) || ~isequal(logical(prepIsGo(:)),isGo(:))), return; end
  rb=S.locoTrace; lens=cellfun(@numel,rb);
  val=lens>=minFrames & ~cellfun(@(v) isempty(v)||all(~isfinite(v)),rb);
  if nnz(val)<6, return; end
  nF=min(lens(val)); gi=find(val(:))'; M=NaN(numel(uTr),nF);
  for a=gi, v=double(rb{a}(:))'; v=v(1:nF); v(~isfinite(v))=0; M(a,:)=v; end
  Mlp=NaN(size(M)); Mlp(gi,:)=lowpass(M(gi,:)',lowpassHz,fs)';
  sd=std(Mlp(gi,:),0,'all','omitnan'); if sd>0, Mn=Mlp/sd; else, Mn=Mlp; end
  if size(Mn,2)>=Ncommon, Mn=Mn(:,1:Ncommon); else, Mn(:,end+1:Ncommon)=NaN; end
  xval=(0:Ncommon-1)*dt-stimOnsetSec; valid=false(numel(uTr),1); valid(gi)=true;
end
