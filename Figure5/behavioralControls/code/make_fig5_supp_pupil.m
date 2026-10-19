%% make_fig5_supp_pupil.m  —  Figures S8-S9: pupil splits on capacity, BOTH windows
% Run for WIN='baseline' [-1.5 0]s and 'post' [1 2]s.
%  - pupil split-capacity (T1/T3/EQ, per WIN) from cache/pupil/
%  - all-trials capacity from the Figure4 capacityTimecourse cache
%  - pupil traces from cache/pupil/<sess>_pupilTraces.mat
%  - splits recomputed deterministically (interp->2Hz LP->window mean->z->
%    within-class tertiles; EQ 1:1 NN Go<->NoGo caliper 0.2 seed 7)
clear; clc; close all;
codeDir=fileparts(mfilename('fullpath')); anaDir=fileparts(codeDir);
pupCache=fullfile(anaDir,'cache','pupil');
allCapDir=fullfile(fileparts(fileparts(anaDir)),'Figure4','manifoldCapacity','cache','capacityTimecourse');
qcDir=fullfile(anaDir,'cache','pupilQC');
panelDir=fullfile(anaDir,'panels','supp'); statsDir=fullfile(anaDir,'stats');
if ~exist(panelDir,'dir'), mkdir(panelDir); end; if ~exist(statsDir,'dir'), mkdir(statsDir); end

DB='BatuhanAngDist'; goVal=45; days={'D1','D2','D3','D4','D5','D6'}; noGoVal=[135,90,75,70,65,60];
dt=0.0646; fs=1/dt; stimOnsetSec=1.5; lowpassHz=2; minFrames=60; evokedWin=[1 2]; xlimSec=[-1 4]; Ncommon=90;
exampleSess='C05_121020';
qc=struct('maxNaNFracWindow',0.10,'jumpMADk',8,'maxJumpFracTrial',0.02,'minStdRaw',1e-6);
eqCaliperSD=0.2; matchSeed=7; strictExclude={'F10_42521_'};
colT1=[0.85 0.15 0.15]; colT3=[0.10 0.30 0.85]; colEQ=[0.50 0.15 0.60]; colAll=[0.15 0.15 0.15];
colGo=[0.20 0.45 0.80]; colNo=[0.85 0.30 0.25]; lsGo='--'; lsNo=':';

QCt=readtable(fullfile(qcDir,sprintf('pupilQC_%s.csv',DB)));
work=QCt(QCt.pass_baseline==1,{'Experiment','BehavCond'});
work=work(~ismember(string(work.Experiment),string(strictExclude)),:);

fid=fopen(fullfile(statsDir,'fig5_supp_pupil_stats.txt'),'w');
for WINc={'baseline','post'}
  WIN=WINc{1};
  if strcmp(WIN,'baseline'), winSec=[-1.5 0]; else, winSec=[1 2]; end
  capT1=cell(6,1);capT3=cell(6,1);PT1Go=cell(6,1);PT1No=cell(6,1);PT3Go=cell(6,1);PT3No=cell(6,1);
  capEQ=cell(6,1);capAll=cell(6,1);EQGo=cell(6,1);EQNo=cell(6,1); tvec=[]; bPup=[];bCap=[];bAllEq=[];
  [exMn,exX,exIsGo,exValid]=computePupilNAS(pupCache,exampleSess,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,Ncommon,winSec,qc);
  exGo=exMn(exValid&exIsGo,:); exNo=exMn(exValid&~exIsGo,:); exGoM=mean(exGo,1,'omitnan'); exNoM=mean(exNo,1,'omitnan');
  for r=1:height(work)
    se=char(work.Experiment(r)); day=char(work.BehavCond(r)); d=find(strcmp(days,day));
    pc=@(part) fullfile(pupCache,sprintf('%s_%s_pupil_%s_%s_capacity_timecourse.mat',se,DB,WIN,part));
    c1=pc('T1'); c3=pc('T3'); ceq=pc('EQ');
    cal=fullfile(allCapDir,sprintf('%s_%s_capacity_timecourse.mat',se,DB));
    haveT1T3=isfile(c1)&&isfile(c3); haveEQ=isfile(ceq)&&isfile(cal);
    [Mn,xval,isGo,valid,metric,usable]=computePupilNAS(pupCache,se,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,Ncommon,winSec,qc);
    if isempty(Mn), continue; end
    zu=metric(usable); z=(metric-mean(zu,'omitnan'))/std(zu,'omitnan');
    [T1,T3]=withinClassTertiles(z,usable,isGo); EQ=matchGoNoGo(z,usable,isGo,eqCaliperSD,matchSeed);
    subT=@(ix) mean(Mn(ix(valid(ix)),:),1,'omitnan');
    if haveT1T3
      C1=load(c1); C3=load(c3); tvec=C1.tvec(:)'; cap1=C1.capacity_rolling(:)'; cap3=C3.capacity_rolling(:)';
      capT1{d}=[capT1{d};cap1]; capT3{d}=[capT3{d};cap3];
      PT1Go{d}=[PT1Go{d};subT(T1(isGo(T1)))]; PT1No{d}=[PT1No{d};subT(T1(~isGo(T1)))];
      PT3Go{d}=[PT3Go{d};subT(T3(isGo(T3)))]; PT3No{d}=[PT3No{d};subT(T3(~isGo(T3)))];
      wl=xval>=winSec(1)&xval<=winSec(2); perTrial=mean(Mn(:,wl),2,'omitnan');
      p1=mean(perTrial(T1(valid(T1))),'omitnan'); p3=mean(perTrial(T3(valid(T3))),'omitnan');
      bPup(end+1,1)=100*(p3-p1)/(0.5*(p3+p1));
      wi=tvec>=evokedWin(1)&tvec<=evokedWin(2); den=0.5*(cap3+cap1); pv=100*(cap3-cap1)./den; pv(~isfinite(pv)|den==0)=NaN; bCap(end+1,1)=mean(pv(wi),'omitnan');
    end
    if haveEQ
      CE=load(ceq); CA=load(cal); if isempty(tvec), tvec=CE.tvec(:)'; end
      bi=tvec<0; ev=tvec>=evokedWin(1)&tvec<=evokedWin(2);
      eqA=CE.capacity_rolling(:)'-mean(CE.capacity_rolling(bi),'omitnan'); alA=CA.capacity_rolling(:)'-mean(CA.capacity_rolling(bi),'omitnan');
      capEQ{d}=[capEQ{d};eqA]; capAll{d}=[capAll{d};alA];
      EQGo{d}=[EQGo{d};subT(EQ(isGo(EQ)))]; EQNo{d}=[EQNo{d};subT(EQ(~isGo(EQ)))];
      evEQ=mean(eqA(ev),'omitnan'); evAL=mean(alA(ev),'omitnan'); bAllEq(end+1,1)=100*(evAL-evEQ)/(0.5*(evAL+evEQ));
    end
  end
  fprintf('WIN=%-8s pupil n=%d med=%.0f%% | cap n=%d med=%.1f%% | ALL-EQ med=%.1f%%\n',WIN,numel(bPup),median(bPup,'omitnan'),numel(bCap),median(bCap,'omitnan'),median(bAllEq(isfinite(bAllEq)),'omitnan'));
  fprintf(fid,'WIN=%s [%g %g]s: pupil n=%d med=%.0f%% (p=%.2g) | cap n=%d med=%.1f%% (p=%.2g) | ALL-EQ med=%.1f%%\n',WIN,winSec(1),winSec(2),numel(bPup),median(bPup,'omitnan'),signrank(bPup),numel(bCap),median(bCap,'omitnan'),signrank(bCap),median(bAllEq(isfinite(bAllEq)),'omitnan'));

  %% figure
  fig=figure('Color','w','Position',[20 20 1040 2300]); TL=tiledlayout(fig,6,1,'TileSpacing','compact','Padding','compact');
  title(TL,sprintf('%s | Pupil is not a confound for discrimination capacity | split=%s [%g %g]s (n=%d)',DB,WIN,winSec(1),winSec(2),numel(bCap)),'Interpreter','none','FontWeight','bold','FontSize',13);
  tTop=tiledlayout(TL,1,2,'TileSpacing','compact','Padding','compact'); tTop.Layout.Tile=1; tTop.Layout.TileSpan=[2 1];
  axA=nexttile(tTop); hold(axA,'on'); plotSingles(axA,exX,exGo,colGo); plotSingles(axA,exX,exNo,colNo);
  plot(axA,exX,exGoM,'-','Color',colGo,'LineWidth',2.6); plot(axA,exX,exNoM,'-','Color',colNo,'LineWidth',2.6);
  xline(axA,0,'k:'); xline(axA,3,'k:'); xlim(axA,xlimSec); xlabel(axA,'Time from stim (s)'); ylabel(axA,'Pupil (AU, /sessStd)');
  title(axA,sprintf('(A) Example %s: single trials + mean (Go blue / NoGo red)',exampleSess),'Interpreter','none'); box(axA,'off'); grid(axA,'on'); axis(axA,'square');
  axD=nexttile(tTop); hold(axD,'on'); cats=categorical({'Pupil T3-T1','Cap T3-T1','Cap ALL-EQ'},{'Pupil T3-T1','Cap T3-T1','Cap ALL-EQ'});
  vAE=bAllEq(isfinite(bAllEq)); X=[repmat(cats(1),numel(bPup),1);repmat(cats(2),numel(bCap),1);repmat(cats(3),numel(vAE),1)]; Y=[bPup;bCap;vAE];
  b=boxchart(axD,X,Y,'BoxWidth',0.5); b.BoxFaceColor=[0.85 0.85 0.85]; b.WhiskerLineColor=[0 0 0]; b.MarkerStyle='none'; rng(1);
  sc=@(xc,v) scatter(axD,repmat(double(xc),numel(v),1)+0.12*(rand(numel(v),1)-0.5),v,26,[0.3 0.3 0.3],'filled','MarkerFaceAlpha',0.6);
  sc(cats(1),bPup); sc(cats(2),bCap); sc(cats(3),vAE); yline(axD,0,'k-'); ylabel(axD,'T3-T1 symmetric %');
  title(axD,sprintf('(D) pupil med %.0f%% (p=%.1g) | cap %.1f%% (p=%.2g) | ALL-EQ %.1f%%',median(bPup,'omitnan'),signrank(bPup),median(bCap,'omitnan'),signrank(bCap),median(vAE,'omitnan')),'Interpreter','tex');
  grid(axD,'on'); box(axD,'off'); axis(axD,'square');
  wrapB=tiledlayout(TL,1,3,'TileSpacing','none','Padding','none'); wrapB.Layout.Tile=3; wrapB.Layout.TileSpan=[2 1];
  tB=tiledlayout(wrapB,2,3,'TileSpacing','compact','Padding','compact'); tB.Layout.Tile=1; tB.Layout.TileSpan=[1 2];   % B panels ~2/3 width (less wide)
  title(tB,'(B) T1 (red, low pupil) vs T3 (blue, high pupil): capacity (left, solid) & pupil (right, Go --/ NoGo :)','FontWeight','bold','FontSize',9);
  for d=1:6
    ax=nexttile(tB); hold(ax,'on'); yyaxis(ax,'left');
    for rr=1:size(capT1{d},1), p=plot(ax,tvec,capT1{d}(rr,:),'-','LineWidth',0.5); p.Color=[colT1 0.18]; p.HandleVisibility='off'; end
    for rr=1:size(capT3{d},1), p=plot(ax,tvec,capT3{d}(rr,:),'-','LineWidth',0.5); p.Color=[colT3 0.18]; p.HandleVisibility='off'; end
    if ~isempty(capT1{d}), plot(ax,tvec,mean(capT1{d},1,'omitnan'),'-','Color',colT1,'LineWidth',2.2); end
    if ~isempty(capT3{d}), plot(ax,tvec,mean(capT3{d},1,'omitnan'),'-','Color',colT3,'LineWidth',2.2); end
    ax.YColor=[0.2 0.2 0.2]; if mod(d,3)==1, ylabel(ax,'Capacity'); end
    yyaxis(ax,'right'); [g1,xg]=avgAl(PT1Go{d},dt,stimOnsetSec); n1=avgAl(PT1No{d},dt,stimOnsetSec); g3=avgAl(PT3Go{d},dt,stimOnsetSec); n3=avgAl(PT3No{d},dt,stimOnsetSec);
    if ~isempty(g1), plot(ax,xg,g1,lsGo,'Color',colT1,'LineWidth',1.4); plot(ax,xg,n1,lsNo,'Color',colT1,'LineWidth',1.4); end
    if ~isempty(g3), plot(ax,xg,g3,lsGo,'Color',colT3,'LineWidth',1.4); plot(ax,xg,n3,lsNo,'Color',colT3,'LineWidth',1.4); end
    ax.YColor=[0.2 0.2 0.2]; if mod(d,3)==0, ylabel(ax,'Pupil (AU)'); end
    xline(ax,0,'k:'); xline(ax,3,'k:'); xlim(ax,xlimSec); grid(ax,'on'); box(ax,'off'); title(ax,sprintf('%s | NoGo=%g',days{d},noGoVal(d))); if d==6, xlabel(ax,'Time from stim (s)'); end
  end
  wrapC=tiledlayout(TL,1,3,'TileSpacing','none','Padding','none'); wrapC.Layout.Tile=5; wrapC.Layout.TileSpan=[2 1];
  tC=tiledlayout(wrapC,2,3,'TileSpacing','compact','Padding','compact'); tC.Layout.Tile=1; tC.Layout.TileSpan=[1 2];   % C panels ~2/3 width (less wide)
  title(tC,'(C) EQ split: baseline-aligned capacity EQ (purple) vs ALL trials (black, left) + EQ pupil (Go --/ NoGo :, right)','FontWeight','bold','FontSize',9);
  for d=1:6
    ax=nexttile(tC); hold(ax,'on'); yyaxis(ax,'left');
    if ~isempty(capEQ{d})
      muE=mean(capEQ{d},1,'omitnan'); seE=std(capEQ{d},0,1,'omitnan')./sqrt(size(capEQ{d},1)); muA=mean(capAll{d},1,'omitnan'); seA=std(capAll{d},0,1,'omitnan')./sqrt(size(capAll{d},1));
      fill(ax,[tvec fliplr(tvec)],[muE+seE fliplr(muE-seE)],colEQ,'EdgeColor','none','FaceAlpha',0.12,'HandleVisibility','off');
      fill(ax,[tvec fliplr(tvec)],[muA+seA fliplr(muA-seA)],colAll,'EdgeColor','none','FaceAlpha',0.12,'HandleVisibility','off');
      plot(ax,tvec,muE,'-','Color',colEQ,'LineWidth',2.2); plot(ax,tvec,muA,'-','Color',colAll,'LineWidth',2.2);
    end
    yline(ax,0,'k-'); ax.YColor=[0.2 0.2 0.2]; if mod(d,3)==1, ylabel(ax,'\Delta Capacity'); end
    yyaxis(ax,'right'); [g,xg]=avgAl(EQGo{d},dt,stimOnsetSec); nn=avgAl(EQNo{d},dt,stimOnsetSec);
    if ~isempty(g), plot(ax,xg,g,lsGo,'Color',colGo,'LineWidth',1.4); plot(ax,xg,nn,lsNo,'Color',colNo,'LineWidth',1.4); end
    ax.YColor=[0.2 0.2 0.2]; if mod(d,3)==0, ylabel(ax,'EQ pupil'); end
    xline(ax,0,'k:'); xline(ax,3,'k:'); xlim(ax,xlimSec); grid(ax,'on'); box(ax,'off'); title(ax,sprintf('%s | NoGo=%g',days{d},noGoVal(d))); if d==6, xlabel(ax,'Time from stim (s)'); end
  end
  exportgraphics(fig,fullfile(panelDir,sprintf('fig5_supp_pupil_%s_%s.png',WIN,DB)),'Resolution',180,'BackgroundColor','white');
  exportgraphics(fig,fullfile(panelDir,sprintf('fig5_supp_pupil_%s_%s.pdf',WIN,DB)),'ContentType','vector','BackgroundColor','white');
  fprintf('saved pupil supp (%s)\n',WIN);
end
fclose(fid); type(fullfile(statsDir,'fig5_supp_pupil_stats.txt'));

%% helpers
function plotSingles(ax,x,M,c), for k=1:size(M,1), y=M(k,:); if all(isnan(y)), continue; end, plot(ax,x,y,'Color',[c 0.05],'LineWidth',0.3,'HandleVisibility','off'); end, end
function [mu,xg]=avgAl(traces,dt,stimOnsetSec), if isempty(traces), mu=[]; xg=[]; return; end, maxF=size(traces,2); mu=mean(traces,1,'omitnan'); xg=(0:maxF-1)*dt-stimOnsetSec; end
function [T1,T3]=withinClassTertiles(z,usable,isGo)
  T1=[];T3=[];
  for cls=[true false]
    idx=find(usable&(isGo==cls)); if numel(idx)<3, continue; end
    vals=z(idx); q33=quantile(vals,1/3); q67=quantile(vals,2/3); T1=[T1; idx(vals<=q33)]; T3=[T3; idx(vals>=q67)];
  end
  T1=sort(T1); T3=sort(T3);
end
function EQ=matchGoNoGo(z,usable,isGo,caliperSD,seed)
  rng(seed); goIdx=find(usable&isGo); noIdx=find(usable&~isGo); EQ=[]; if isempty(goIdx)||isempty(noIdx), return; end
  zg=z(goIdx); zn=z(noIdx); usedNo=false(numel(noIdx),1); order=randperm(numel(goIdx));
  for ii=order, dd=abs(zn-zg(ii)); dd(usedNo)=Inf; [dm,kk]=min(dd); if isfinite(dm)&&dm<=caliperSD, usedNo(kk)=true; EQ=[EQ; goIdx(ii); noIdx(kk)]; end, end
  EQ=sort(EQ);
end
function [Mn,xval,isGo,valid,metric,usable]=computePupilNAS(pupDir,se,DB,goVal,dt,fs,stimOnsetSec,lowpassHz,minFrames,Ncommon,winSec,qc)
  Mn=[];xval=[];isGo=[];valid=[];metric=[];usable=[];
  f=fullfile(pupDir,sprintf('%s_%s_pupilTraces.mat',se,DB)); if ~isfile(f), return; end
  S=load(f); uTr=double(S.Trial); isGo=double(S.Visual_Stim)==goVal; rb=S.pupilTrace; n=numel(uTr);
  LP=cell(n,1); metric=NaN(n,1); usable=false(n,1); lens=zeros(n,1);
  for a=1:n
    v=rb{a}; if isempty(v)||all(isnan(v)), continue; end
    v=double(v(:)'); L=numel(v); lens(a)=L; if L<minFrames, continue; end
    tRel=(0:L-1)*dt-stimOnsetSec; idx=tRel>=winSec(1)&tRel<=winSec(2);
    vi=v; if any(isnan(vi)), g=find(~isnan(vi)); if numel(g)<3, continue; end, vi=interp1(g,vi(g),1:L,'linear','extrap'); end
    Wn=lowpassHz/(fs/2); vLP=vi; if Wn<1, [bb,aa]=butter(2,Wn,'low'); vLP=filtfilt(bb,aa,vi); end
    LP{a}=vLP;
    if any(idx)
      cov=mean(~isnan(v(idx))); dd=diff(vi); okjump=true;
      if numel(dd)>=5, m=mad(dd,1); jf=(m>0)*mean(abs(dd-median(dd))>qc.jumpMADk*m); okjump=jf<=qc.maxJumpFracTrial; end
      if cov>=(1-qc.maxNaNFracWindow) && std(vi(idx),0,'omitnan')>qc.minStdRaw && okjump
        metric(a)=mean(vLP(idx),'omitnan'); usable(a)=true;
      end
    end
  end
  gi=find(cellfun(@(x)~isempty(x),LP))'; if numel(gi)<6, return; end
  nF=min(lens(gi)); M=NaN(n,nF); for a=gi, M(a,:)=LP{a}(1:nF); end
  sd=std(M(gi,:),0,'all','omitnan'); if sd>0, M=M/sd; end
  if size(M,2)>=Ncommon, Mn=M(:,1:Ncommon); else, Mn=[M NaN(n,Ncommon-size(M,2))]; end
  xval=(0:Ncommon-1)*dt-stimOnsetSec; valid=false(n,1); valid(gi)=true;
end
