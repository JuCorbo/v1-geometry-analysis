%% make_fig5_supp_supermouse.m  —  Figure S7: super-mouse capacity + lick (All/Lick/No-lick)
%  - capacity + aux from the FIXEDTRIALS summaries (cache/supermouse)
%  - per-session lick grids rebuilt from the cached lick export (lickTracesFull)
%  - pseudo-lick reconstructed from the stored aux draws
% Layout: 6 days, each with capacity (top) + reconstructed lick (bottom).
clear; clc;
codeDir=fileparts(mfilename('fullpath')); anaDir=fileparts(codeDir); figDir=fileparts(fileparts(anaDir));
sumDir  = fullfile(anaDir,'cache','supermouse');
lickDir = fullfile(anaDir,'cache','lickTracesFull');
sharedTr= fullfile(figDir,'_sharedWindowedData','data','trained');
panelDir= fullfile(anaDir,'panels','supp');
DB='BatuhanAngDist'; days={'D1','D2','D3','D4','D5','D6'}; nogo=[135 90 75 70 65 60];
stimOnsetSec=1.5; doSmooth=true; smoothWinSec=0.15; lickPlotEndSec=4;
colAll=[0 0 0]; colLick=[0 0.6 0]; colNo=[0.85 0 0];

% day -> alphabetically-sorted session list (matches the compute's dir() order)
daySess=cell(6,1);
for d=1:6
  ff=dir(fullfile(sharedTr,sprintf('D%d',d),'*_win.mat'));
  s=erase(erase({ff.name},[DB '_']),'_win.mat'); daySess{d}=sort(s);
end

fig=figure('Color','w','Position',[40 40 1500 950]);
tl=tiledlayout(fig,4,3,'Padding','compact','TileSpacing','compact');
title(tl,'Super-mouse: capacity (rows 1,3) + reconstructed lick (rows 2,4) | All (black) / Lick (green) / No-lick (red)','FontWeight','bold');
order=[1 2 3 1 2 3 4 5 6 4 5 6];  % day per tile (cap row then lick row per block)
for blk=0:1
  for c=1:3
    d=blk*3+c;
    S=load(fullfile(sumDir,sprintf('supermouse_%s_FIXEDTRIALS_summary.mat',days{d})));
    tCap=S.tvec(:)';
    % capacity tile
    axc=nexttile(tl, blk*6 + c); hold(axc,'on');
    band(axc,tCap,S.capAll,colAll); band(axc,tCap,S.capLick,colLick); band(axc,tCap,S.capNoLick,colNo);
    xline(axc,0,'k--','HandleVisibility','off'); xline(axc,3,'k:','HandleVisibility','off');
    xlim(axc,[-1 4]); ylabel(axc,'Capacity'); title(axc,sprintf('%s | NoGo=%d',days{d},nogo(d))); box(axc,'off');
    % lick tile
    [grids,tL]=buildGrids(lickDir,DB,daySess{d},stimOnsetSec,doSmooth,smoothWinSec,lickPlotEndSec);
    lAll=reconAux(S.auxAll,grids); lLick=reconAux(S.auxLick,grids); lNo=reconAux(S.auxNoLick,grids);
    axl=nexttile(tl, blk*6 + 3 + c); hold(axl,'on');
    band(axl,tL,lAll,colAll); band(axl,tL,lLick,colLick); band(axl,tL,lNo,colNo);
    xline(axl,0,'k--','HandleVisibility','off'); xline(axl,3,'k:','HandleVisibility','off');
    xlim(axl,[-1 4]); ylabel(axl,'Lick (a.u.)'); xlabel(axl,'Time from stim (s)'); box(axl,'off');
    if d==1, legend(axc,{'All','Lick','No-lick'},'Location','best','Box','off'); end
  end
end
base=fullfile(panelDir,'fig5_supp_supermouse_capLick');
exportgraphics(fig,[base '.png'],'Resolution',180,'BackgroundColor','white');
exportgraphics(fig,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
fprintf('saved %s\n',base);

%% ---- helpers ----
function h=band(ax,t,M,c)
  mu=mean(M,1,'omitnan'); se=std(M,0,1,'omitnan')./sqrt(sum(~isnan(M),1));
  fill(ax,[t fliplr(t)],[mu+se fliplr(mu-se)],c,'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
  h=plot(ax,t,mu,'-','Color',c,'LineWidth',2.0);
end
function [grids,tCommon]=buildGrids(lickDir,DB,expList,stimOnsetSec,doSmooth,smoothWinSec,lickPlotEndSec)
  nSess=numel(expList); dtSess=nan(nSess,1); decoded=cell(nSess,1);
  for s=1:nSess
    f=fullfile(lickDir,sprintf('%s_%s_lickTraces.mat',expList{s},DB));
    if ~isfile(f), decoded{s}={}; continue; end
    S=load(f); [~,ord]=sort(double(S.Trial));
    vs=double(S.Visual_Stim(ord)); rs=double(S.Response(ord)); dtv=double(S.Lick_dt(ord)); tr=S.lickTrace(ord);
    ok=isfinite(vs)&isfinite(rs); tr=tr(ok); dtv=dtv(ok);
    decoded{s}=tr; d=dtv(isfinite(dtv)&dtv>0); if ~isempty(d), dtSess(s)=median(d); end
  end
  dtCommon=median(dtSess(isfinite(dtSess)));
  maxDur=0; for s=1:nSess, for i=1:numel(decoded{s}), maxDur=max(maxDur,(numel(decoded{s}{i})-1)*dtSess(s)); end, end
  if ~(dtCommon>0)||maxDur<=0, grids=cell(nSess,1); tCommon=[]; return; end
  Lc=floor(maxDur/dtCommon)+1; tAbs=(0:Lc-1)*dtCommon; tCommon=tAbs-stimOnsetSec;
  keep=tCommon<=lickPlotEndSec; tAbs=tAbs(keep); tCommon=tCommon(keep); Lc=numel(tCommon);
  grids=cell(nSess,1);
  for s=1:nSess
    tr=decoded{s}; dts=dtSess(s); if isempty(tr)||~(dts>0), grids{s}=[]; continue; end
    M=nan(numel(tr),Lc);
    for i=1:numel(tr)
      v=double(tr{i}(:))'; if isempty(v), continue; end
      if doSmooth&&smoothWinSec>0, v=movmean(v,max(1,round(smoothWinSec/dts)),'Endpoints','shrink'); end
      t0=(0:numel(v)-1)*dts;
      if abs(dts-dtCommon)>1e-9, v=interp1(t0,v,tAbs,'linear','extrap'); else, if numel(v)<numel(tAbs), v(numel(tAbs))=NaN; end; v=v(1:numel(tAbs)); end
      M(i,:)=v(1:Lc);
    end
    grids{s}=M;
  end
end
function lickBoth=reconAux(aux,grids)
  nRep=numel(aux); L=0;
  for s=1:numel(grids), if ~isempty(grids{s}), L=size(grids{s},2); break; end, end
  lickBoth=nan(nRep,L);
  for r=1:nRep
    A=aux{r}; if ~isfield(A,'draw')||~isfield(A,'nNeurBySess')||~isfield(A.draw,'sessIdxVec'), continue; end
    draw=A.draw; sIdx=draw.sessIdxVec(:)'; w=A.nNeurBySess(:)';
    Wtot=sum(w(sIdx)); if ~(Wtot>0), Wtot=1; end
    nP=size(draw.trialsGoUsed,1); accG=zeros(nP,L); accN=zeros(nP,L);
    for p=1:nP
      aG=zeros(1,L); aN=zeros(1,L);
      for ii=1:numel(sIdx)
        s=sIdx(ii); M=grids{s}; if isempty(M), continue; end
        tiG=draw.trialsGoUsed(p,ii); tiN=draw.trialsNoGoUsed(p,ii);
        if isfinite(tiG)&&tiG>=1&&tiG<=size(M,1), v=M(tiG,:); v(~isfinite(v))=0; aG=aG+w(s)*v; end
        if isfinite(tiN)&&tiN>=1&&tiN<=size(M,1), v=M(tiN,:); v(~isfinite(v))=0; aN=aN+w(s)*v; end
      end
      accG(p,:)=aG/Wtot; accN(p,:)=aN/Wtot;
    end
    lickBoth(r,:)=0.5*(mean(accG,1,'omitnan')+mean(accN,1,'omitnan'));
  end
end
