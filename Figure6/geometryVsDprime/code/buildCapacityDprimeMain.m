%% buildCapacityDprimeMain.m  (Figure 6 — panels c, d, e: capacity vs d')
%    c: correlation timecourse rho(t) of capacity vs d', raw (solid) and with
%       recording day/angle partialled out (dashed)
%    d: 1-2 s scatter, full (capacity vs d', all sessions; pooled Spearman)
%    e: 1-2 s scatter with difficulty regressed out (within-day-centred capacity
%       vs within-day-centred d'; within-angle partial Spearman)
%  SQL-free: d' from behavDprimeBySession (deposit CSV); capacity from local caches.
clear; clc;
goStim=45; scatterWin=[1 2];
thisDir=fileparts(mfilename('fullpath')); if isempty(thisDir),thisDir=pwd; end
anaDir=fileparts(thisDir); finalRoot=fileparts(fileparts(anaDir));
addpath(fullfile(finalRoot,'_sharedBehavior','code'));
capDir=fullfile(finalRoot,'Figure4','manifoldCapacity','cache','capacityTimecourse');
mainDir=fullfile(anaDir,'panels','main'); if ~exist(mainDir,'dir'),mkdir(mainDir); end
dataDir=fullfile(anaDir,'data'); if ~exist(dataDir,'dir'),mkdir(dataDir); end
statsDir=fullfile(anaDir,'stats'); if ~exist(statsDir,'dir'),mkdir(statsDir); end

[dprime,~,~,sessAllDays]=behavDprimeBySession('Trained',goStim,10);
capMap=containers.Map;
dd=dir(fullfile(capDir,'*_BatuhanAngDist_capacity_timecourse.mat'));
for f=1:numel(dd), id=extractBefore(dd(f).name,'_BatuhanAngDist_capacity_timecourse.mat');
  if ~isempty(id), capMap(id)=load(fullfile(dd(f).folder,dd(f).name),'capacity_rolling','tvec'); end; end

tCommon=[]; Cmat=[]; Dvec=[]; DayVec=[];
for day=1:6, ids=sessAllDays{day};
  for a=1:numel(ids), id=ids{a}; if ~isKey(capMap,id), continue; end
    S=capMap(id); t=S.tvec(:); c=S.capacity_rolling(:);
    if isempty(tCommon), tCommon=t; end
    if numel(t)~=numel(tCommon)||any(abs(t-tCommon)>1e-9), c=interp1(t,c,tCommon,'linear','extrap'); end
    Cmat(end+1,:)=c'; Dvec(end+1,1)=dprime(day,a); DayVec(end+1,1)=day; %#ok<AGROW>
  end
end
Ddum=zeros(numel(DayVec),5); for d=2:6, Ddum(:,d-1)=(DayVec==d); end
dayCol=lines(6);

% raw correlation timecourse
nT=numel(tCommon); rhoRaw=nan(1,nT);
for ti=1:nT, x=Cmat(:,ti); ok=isfinite(x)&isfinite(Dvec);
  if sum(ok)>=8, rhoRaw(ti)=corr(x(ok),Dvec(ok),'type','Spearman'); end
end

% 1-2 s window capacity
widx=tCommon>=scatterWin(1)&tCommon<scatterWin(2); xW=mean(Cmat(:,widx),2,'omitnan');
ok=isfinite(xW)&isfinite(Dvec);
[rFull,pFull]=corr(xW(ok),Dvec(ok),'type','Spearman');
[rReg,pReg]=partialcorr(xW(ok),Dvec(ok),Ddum(ok,:),'type','Spearman');
xc=xW; yc=Dvec; for d=1:6, m=DayVec==d; xc(m)=xW(m)-mean(xW(m),'omitnan'); yc(m)=Dvec(m)-mean(Dvec(m),'omitnan'); end

fig=figure('Color','w','Position',[80 100 1350 400]); tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% (a) timecourse
ax=nexttile(tl,1); hold(ax,'on'); box(ax,'off');
patch(ax,[scatterWin(1) scatterWin(2) scatterWin(2) scatterWin(1)],[-1 -1 1 1],[0.9 0.9 0.6],'EdgeColor','none','FaceAlpha',0.4);
plot(ax,tCommon,rhoRaw,'k-','LineWidth',2);
yline(ax,0,':k'); xline(ax,0,'-','Color',[0.5 0.5 0.5]);
xlabel(ax,'Time from stimulus onset (s)'); ylabel(ax,'\rho (capacity vs d'')');
ylim(ax,[-0.2 0.7]);
title(ax,'a  Correlation timecourse','Interpreter','none');

% (b) full 1-2 s scatter
ax=nexttile(tl,2); hold(ax,'on'); box(ax,'off'); axis(ax,'square');
for d=1:6, m=ok&DayVec==d; scatter(ax,xW(m),Dvec(m),44,dayCol(d,:),'filled','MarkerFaceAlpha',0.8); end
pf=polyfit(xW(ok),Dvec(ok),1); xr=linspace(min(xW(ok)),max(xW(ok)),50); plot(ax,xr,polyval(pf,xr),'k-','LineWidth',2);
xlabel(ax,'Capacity (1-2 s)'); ylabel(ax,'d''');
title(ax,sprintf('b  Full: \\rho=%.2f, p=%.1e',rFull,pFull),'Interpreter','tex');

% (c) difficulty-regressed 1-2 s scatter
ax=nexttile(tl,3); hold(ax,'on'); box(ax,'off'); axis(ax,'square');
for d=1:6, m=ok&DayVec==d; scatter(ax,xc(m),yc(m),44,dayCol(d,:),'filled','MarkerFaceAlpha',0.8); end
pf=polyfit(xc(ok),yc(ok),1); xr=linspace(min(xc(ok)),max(xc(ok)),50); plot(ax,xr,polyval(pf,xr),'k-','LineWidth',2);
xline(ax,0,':','Color',0.7*[1 1 1]); yline(ax,0,':','Color',0.7*[1 1 1]);
xlabel(ax,'Capacity (1-2 s), within-angle'); ylabel(ax,'d'', within-angle');
title(ax,sprintf('c  Difficulty regressed: \\rho=%.2f, p=%.2g',rReg,pReg),'Interpreter','tex');

exportgraphics(fig,fullfile(mainDir,'fig6_cde_capacityVsDprime.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(mainDir,'fig6_cde_capacityVsDprime.png'),'Resolution',200);

fid=fopen(fullfile(statsDir,'capacityVsDprime_stats.txt'),'w');
fprintf(fid,'Capacity vs d'' (main). n=%d trained sessions.\n',sum(ok));
fprintf(fid,'1-2 s window:  FULL Spearman rho=%.3f p=%.3g | DIFFICULTY-REGRESSED (partial | day) rho=%.3f p=%.3g\n',rFull,pFull,rReg,pReg);
fprintf(fid,'raw association is carried by the between-day (Go/NoGo angle) axis; within angle it is null.\n');
fclose(fid);
fprintf('FULL rho=%.2f (p=%.2g) | REGRESSED rho=%.2f (p=%.2g)\n',rFull,pFull,rReg,pReg);
fprintf('saved fig6_capacityVsDprime_main\n');
