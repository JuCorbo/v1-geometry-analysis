%% buildPacmapSepVsDprime.m  (Figure S12 — PaCMAP Go/NoGo separation vs d')
%  Separation metric: spread-normalised Go-to-NoGo centroid distance in the main
%  global-z PaCMAP embedding (seed0):
%       centDistNorm = ||cGo - cNo|| / mean(rGo, rNo)
%  (rGo, rNo = mean distance of single trials to their own class centroid; the
%  normalisation removes the arbitrary per-session PaCMAP scale). All windows of
%  a session share one embedding frame, so the metric is averaged over the
%  1-2 s stimulus windows.
%
%  Panels:
%    a: correlation timecourse rho(t) of centDistNorm vs d'
%    b: 1-2 s scatter, full (pooled Spearman)
%    c: 1-2 s scatter with difficulty regressed out (within-day-centred;
%       partial Spearman controlling for recording day)
%  SQL-free: d' from behavDprimeBySession (deposit CSV); geometry from the Fig 2
%  PaCMAP caches. A per-session metric table is cached to ../data for re-plotting.
clear; clc;
goStim=45; scatterWin=[1 2]; seed=0; dims=[1 2 3];
thisDir=fileparts(mfilename('fullpath')); if isempty(thisDir), thisDir=pwd; end
anaDir=fileparts(thisDir); finalRoot=fileparts(fileparts(anaDir));
addpath(fullfile(finalRoot,'_sharedBehavior','code'));
embRoot=fullfile(finalRoot,'Figure2','pacmap','cache','main_pacmap_bwgz_pm080', ...
                 sprintf('baseline_seed%d',seed),'BatuhanAngDist');
assert(isfolder(embRoot),'Embedding root not found:\n  %s',embRoot);
suppDir=fullfile(anaDir,'panels','supp'); if ~exist(suppDir,'dir'), mkdir(suppDir); end
dataDir=fullfile(anaDir,'data');           if ~exist(dataDir,'dir'), mkdir(dataDir); end
statsDir=fullfile(anaDir,'stats');         if ~exist(statsDir,'dir'), mkdir(statsDir); end

% d' per session -> map experimentID -> d'
[dprime,~,~,sessAllDays]=behavDprimeBySession('Trained',goStim,10);
dMap=containers.Map;
for day=1:6, ids=sessAllDays{day};
  for a=1:numel(ids), if ~isempty(ids{a}), dMap(ids{a})=dprime(day,a); end; end
end

files=[]; for d=1:6, files=[files; dir(fullfile(embRoot,sprintf('D%d',d),'*_PaCMAP_UMATO.mat'))]; end %#ok<AGROW>
assert(~isempty(files),'No embedding files under %s',embRoot);

tCommon=[]; rows=struct('exp',{},'day',{},'dp',{},'cdn',{},'cdnT',{});
for fi=1:numel(files)
  S=load(fullfile(files(fi).folder,files(fi).name));
  if ~all(isfield(S,{'Y_pacmap','labels','trial_ids','time_window_idx','time_window_center'})), continue; end
  Y=double(S.Y_pacmap(:,dims)); lab=double(S.labels(:));
  twi=double(S.time_window_idx(:)); twc=double(S.time_window_center(:));
  expID=getExperimentID(S,files(fi).name);
  day=str2double(regexp(files(fi).folder,'D(\d)$','tokens','once'));
  if ~isKey(dMap,expID), continue; end
  sv=unique(lab); sv=sv(isfinite(sv)); if numel(sv)~=2||~any(sv==goStim), continue; end
  nogo=sv(sv~=goStim); nogo=nogo(1);
  uW=sort(unique(twi(isfinite(twi)))); tC=arrayfun(@(w)mean(twc(twi==w),'omitnan'),uW);
  if isempty(tCommon), tCommon=tC(:)'; end
  nW=numel(uW); cdnT=nan(1,nW);
  for wi=1:nW
    m=(twi==uW(wi)); Xg=Y(m&lab==goStim,:); Xn=Y(m&lab==nogo,:);
    if size(Xg,1)<3||size(Xn,1)<3, continue; end
    cg=mean(Xg,1); cn=mean(Xn,1);
    rg=mean(sqrt(sum((Xg-cg).^2,2))); rn=mean(sqrt(sum((Xn-cn).^2,2)));
    cdnT(wi)=norm(cg-cn)/max(0.5*(rg+rn),eps);
  end
  widx=tC>=scatterWin(1)&tC<scatterWin(2);
  rows(end+1)=struct('exp',expID,'day',day,'dp',dMap(expID), ...
    'cdn',mean(cdnT(widx),'omitnan'),'cdnT',cdnT); %#ok<AGROW>
end
fprintf('sessions used: %d\n',numel(rows));

Dvec=[rows.dp]'; DayVec=[rows.day]'; xW=[rows.cdn]';
Ddum=zeros(numel(DayVec),5); for d=2:6, Ddum(:,d-1)=(DayVec==d); end
dayCol=lines(6);

% raw correlation timecourse
Xt=cat(1,rows.cdnT); nT=numel(tCommon); rhoRaw=nan(1,nT);
for ti=1:nT, x=Xt(:,ti); ok=isfinite(x)&isfinite(Dvec);
  if sum(ok)>=8, rhoRaw(ti)=corr(x(ok),Dvec(ok),'type','Spearman'); end
end

% 1-2 s window stats
ok=isfinite(xW)&isfinite(Dvec);
[rFull,pFull]=corr(xW(ok),Dvec(ok),'type','Spearman');
[rReg,pReg]  =partialcorr(xW(ok),Dvec(ok),Ddum(ok,:),'type','Spearman');
% within-day-centred values for the difficulty-regressed scatter
xc=xW; yc=Dvec;
for d=1:6, m=DayVec==d; xc(m)=xW(m)-mean(xW(m),'omitnan'); yc(m)=Dvec(m)-mean(Dvec(m),'omitnan'); end

fig=figure('Color','w','Position',[80 100 1350 400]);
tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% (a) raw timecourse
ax=nexttile(tl,1); hold(ax,'on'); box(ax,'off');
patch(ax,[scatterWin(1) scatterWin(2) scatterWin(2) scatterWin(1)],[-1 -1 1 1],[0.9 0.9 0.6],'EdgeColor','none','FaceAlpha',0.4);
plot(ax,tCommon,rhoRaw,'k-','LineWidth',2);
yline(ax,0,':k'); xline(ax,0,'-','Color',[0.5 0.5 0.5]);
xlabel(ax,'Time from stimulus onset (s)'); ylabel(ax,'\rho (separation vs d'')');
ylim(ax,[-0.3 0.8]); title(ax,'a  Correlation timecourse','Interpreter','none');

% (b) full 1-2 s scatter
ax=nexttile(tl,2); hold(ax,'on'); box(ax,'off'); axis(ax,'square');
for d=1:6, m=ok&DayVec==d; scatter(ax,xW(m),Dvec(m),44,dayCol(d,:),'filled','MarkerFaceAlpha',0.8); end
pf=polyfit(xW(ok),Dvec(ok),1); xr=linspace(min(xW(ok)),max(xW(ok)),50); plot(ax,xr,polyval(pf,xr),'k-','LineWidth',2);
xlabel(ax,'Normalised centroid distance (1-2 s)'); ylabel(ax,'d''');
title(ax,sprintf('b  Full: \\rho=%.2f, p=%.1e',rFull,pFull),'Interpreter','tex');

% (c) difficulty-regressed 1-2 s scatter
ax=nexttile(tl,3); hold(ax,'on'); box(ax,'off'); axis(ax,'square');
for d=1:6, m=ok&DayVec==d; scatter(ax,xc(m),yc(m),44,dayCol(d,:),'filled','MarkerFaceAlpha',0.8); end
pf=polyfit(xc(ok),yc(ok),1); xr=linspace(min(xc(ok)),max(xc(ok)),50); plot(ax,xr,polyval(pf,xr),'k-','LineWidth',2);
xline(ax,0,':','Color',0.7*[1 1 1]); yline(ax,0,':','Color',0.7*[1 1 1]);
xlabel(ax,'Norm. centroid distance (1-2 s), within-angle'); ylabel(ax,'d'', within-angle');
title(ax,sprintf('c  Difficulty regressed: \\rho=%.2f, p=%.2g',rReg,pReg),'Interpreter','tex');

exportgraphics(fig,fullfile(suppDir,'pacmapSepVsDprime.pdf'),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(suppDir,'pacmapSepVsDprime.png'),'Resolution',200);

% cache per-session table + stats
Tsess=table(string({rows.exp})',[rows.day]',[rows.dp]',[rows.cdn]', ...
  'VariableNames',{'Experiment','Day','dprime','centDistNorm_1to2s'});
save(fullfile(dataDir,'pacmapSepVsDprime_sessions.mat'),'Tsess','tCommon','rhoRaw','seed','scatterWin');
fid=fopen(fullfile(statsDir,'pacmapSepVsDprime_stats.txt'),'w');
fprintf(fid,'PaCMAP Go/NoGo separation vs d'' (spread-normalised centroid distance). n=%d trained sessions.\n',sum(ok));
fprintf(fid,'metric = ||cGo-cNo|| / mean(rGo,rNo), main global-z PaCMAP seed0, averaged over 1-2 s windows.\n');
fprintf(fid,'1-2 s window:  FULL Spearman rho=%.3f p=%.3g | DIFFICULTY-REGRESSED (partial | day) rho=%.3f p=%.3g\n',rFull,pFull,rReg,pReg);
fprintf(fid,'the association survives day-partialling (unlike manifold capacity), i.e. it retains a within-difficulty link to d''.\n');
fclose(fid);
fprintf('FULL rho=%.2f (p=%.2g) | REGRESSED rho=%.2f (p=%.2g)\n',rFull,pFull,rReg,pReg);
fprintf('saved fig6_pacmapSepVsDprime_main\n');

%% helpers
function experimentID=getExperimentID(S,fname)
  experimentID='';
  if isfield(S,'experimentID')&&~isempty(S.experimentID), experimentID=matToChar(S.experimentID); end
  if isempty(experimentID)
    [~,base]=fileparts(fname); base=regexprep(base,'_PaCMAP_UMATO$',''); base=regexprep(base,'^BatuhanAngDist_','');
    experimentID=base;
  end
end
function out=matToChar(x)
  while iscell(x)&&numel(x)==1, x=x{1}; end
  if isnumeric(x), out=strtrim(char(x(:)')); return; end
  if isstring(x), x=char(x); end
  if ischar(x)
    t=strtrim(x); nums=regexp(t,'-?\d+','match');
    if ~isempty(nums), v=str2double(nums);
      if numel(v)>1&&all(isfinite(v))&&all((v>=32&v<=126)|v==9|v==10|v==13), out=strtrim(char(v(:)')); return; end
    end
    out=t; return;
  end
  out=strtrim(char(string(x)));
end
