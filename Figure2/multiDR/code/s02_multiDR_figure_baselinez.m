%% multiDRfigure.m  (Stage B: two-pass align + discriminant-view cloud figure)
% Loads cached embeddings, aligns each method/group EXACTLY like
% pacmapThroughDaysPooled.m (within-day medoid -> day->global, Go+NoGo
% centroid-trajectory landmarks), then displays clouds in a 2D DISCRIMINANT view
% that maximizes Go/NoGo separation: x = Go->NoGo separation axis, y = top
% orthogonal variance direction. The SAME projection is used across all 6 panels
% of a method, so trained-cloud consistency across days is directly visible.
% Zoomed to the cloud bulk (2-98 pct).

clear; clc;
codeDir=fileparts(mfilename('fullpath'));
cacheFile=getenv('MULTIDR_CACHE');
if isempty(cacheFile), cacheFile=fullfile(fileparts(codeDir),'cache','multiDR_embeddings_baselinez.mat'); end
L=load(cacheFile); E=L.E;
snapT=2.5; daysToPlot=1:6;
methods={'Ypac','PaCMAP';'Yuma','UMATO';'Ycmds','classical MDS';'Yiso','Isomap';'Ynmds','nMDS'};
colGo=[0.20 0.40 0.80]; colNo=[0.80 0.20 0.20]; mkSz=22;
outDir=getenv('MULTIDR_OUT');
if isempty(outDir), outDir=fullfile(fileparts(codeDir),'panels','supp'); end
if ~exist(outDir,'dir'), mkdir(outDir); end
grp=string({E.group})'; nD=numel(daysToPlot); cols=struct('Naive',1:nD,'Trained',nD+1:2*nD);

% --- Phase 1: align, discriminant projection (per method); zoom limits PER (method,group) ---
PD=cell(size(methods,1),1);   % PD{mi}.(group){dayidx}={qg,qn}; PD{mi}.lim.(group)=[xl;yl]
for mi=1:size(methods,1)
  fld=methods{mi,1}; GA=struct(); PgAll=[]; PnAll=[];
  for gg=["Naive","Trained"]
    S=struct('Y',{},'labels',{},'twi',{},'twc',{},'goV',{},'noV',{},'day',{});
    for j=find(grp==gg)'
      Y=E(j).(fld); if isempty(Y), continue; end
      S(end+1)=struct('Y',Y,'labels',E(j).labels,'twi',E(j).twi,'twc',E(j).twc,'goV',E(j).goV,'noV',E(j).noV,'day',E(j).day); %#ok<SAGROW>
    end
    A=twoPassAlignGroup(S,snapT); GA.(gg)=A;
    for a=1:numel(A)
      PgAll=[PgAll; A(a).cloudY(A(a).cloudLab==A(a).goV,:)]; %#ok<AGROW>
      PnAll=[PnAll; A(a).cloudY(A(a).cloudLab==A(a).noV,:)]; %#ok<AGROW>
    end
  end
  if isempty(PgAll)||isempty(PnAll), PD{mi}=[]; continue; end
  u1=(mean(PnAll,1)-mean(PgAll,1))'; if norm(u1)<eps, u1=[1;zeros(size(PgAll,2)-1,1)]; end; u1=u1/norm(u1);
  allP=[PgAll;PnAll]; res=allP-(allP*u1)*u1'; cf=pca(res);
  if isempty(cf), u2=null(u1'); u2=u2(:,1); else, u2=cf(:,1); end
  W=[u1 u2]; proj=@(P)P*W;
  pan=struct(); lim=struct();
  for gg=["Naive","Trained"]
    A=GA.(gg); gAcc=[];
    for ci=1:numel(daysToPlot)
      d=daysToPlot(ci); Pg=[]; Pn=[];
      for a=1:numel(A)
        if A(a).day~=d, continue; end
        Pg=[Pg; A(a).cloudY(A(a).cloudLab==A(a).goV,:)]; Pn=[Pn; A(a).cloudY(A(a).cloudLab==A(a).noV,:)]; %#ok<AGROW>
      end
      qg=[]; qn=[]; if ~isempty(Pg), qg=proj(Pg); end; if ~isempty(Pn), qn=proj(Pn); end
      pan.(gg){ci}={qg,qn}; gAcc=[gAcc; qg; qn]; %#ok<AGROW>
    end
    xl=prctile(gAcc(:,1),[1.5 98.5]); yl=prctile(gAcc(:,2),[1.5 98.5]);
    xl=xl+0.12*[-1 1]*diff(xl); yl=yl+0.12*[-1 1]*diff(yl);
    lim.(gg)=[xl; yl];
  end
  pan.lim=lim; PD{mi}=pan;
end

% --- Phase 2: plot (zoom shared within each method+group only) ---
fig=figure('Color','w','Position',[20 20 2100 1450]);
tl=tiledlayout(size(methods,1),2*nD,'TileSpacing','compact','Padding','compact');
for mi=1:size(methods,1)
  pan=PD{mi}; if isempty(pan), continue; end
  for gg=["Naive","Trained"]
    lg=pan.lim.(gg);
    for ci=1:numel(daysToPlot)
      d=daysToPlot(ci); ax=nexttile(tl,(mi-1)*2*nD+cols.(gg)(ci)); hold(ax,'on'); box(ax,'on');
      q=pan.(gg){ci}; qg=q{1}; qn=q{2};
      if ~isempty(qg), scatter(ax,qg(:,1),qg(:,2),mkSz,colGo,'filled','MarkerFaceAlpha',0.07,'MarkerEdgeColor','none'); end
      if ~isempty(qn), scatter(ax,qn(:,1),qn(:,2),mkSz,colNo,'filled','MarkerFaceAlpha',0.07,'MarkerEdgeColor','none'); end
      if size(qg,1)>=5, e=ellipseXY(qg,1); plot(ax,e(1,:),e(2,:),'-','Color',colGo,'LineWidth',2); end
      if size(qn,1)>=5, e=ellipseXY(qn,1); plot(ax,e(1,:),e(2,:),'-','Color',colNo,'LineWidth',2); end
      xlim(ax,lg(1,:)); ylim(ax,lg(2,:)); set(ax,'XTickLabel',[],'YTickLabel',[]); grid(ax,'on');
      if mi==1, title(ax,sprintf('%s D%d',gg,d),'FontSize',10); end
      if cols.(gg)(ci)==1, ylabel(ax,methods{mi,2},'FontWeight','bold','FontSize',11); end
    end
  end
end
sgtitle(tl,'Discriminant view: faint clouds + 1-SD covariance ellipses (Go blue, NoGo red); zoom shared within method+group');
exportgraphics(fig,fullfile(outDir,'multiDR_clouds_ellipse.png'),'Resolution',150,'BackgroundColor','white');
try, exportgraphics(fig,fullfile(outDir,'multiDR_clouds_ellipse.pdf'),'ContentType','image','BackgroundColor','white'); catch, end
try, savefig(fig,fullfile(outDir,'multiDR_clouds_ellipse.fig')); catch, end
fprintf('Saved multiDR_clouds_ellipse.png to %s\n',outDir);

%% --- Figure 2: per method, Naive & Trained panels with all 3 days' (D2-4) ellipses overlaid ---
blend=@(c,f)(1-f)*c+f*[1 1 1];
shGo=cell(1,nD); shNo=cell(1,nD);                    % light->dark = D1->D6
for d=1:nD, fr=0.62*(nD-d)/max(nD-1,1); shGo{d}=blend(colGo,fr); shNo{d}=blend(colNo,fr); end
fig2=figure('Color','w','Position',[60 20 720 1450]);
tl2=tiledlayout(size(methods,1),2,'TileSpacing','compact','Padding','compact');
groups=["Naive","Trained"];
for mi=1:size(methods,1)
  pan=PD{mi}; if isempty(pan), continue; end
  for gc=1:2
    gg=groups(gc); lg=pan.lim.(gg);
    ax=nexttile(tl2,(mi-1)*2+gc); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
    for ci=1:numel(daysToPlot)
      q=pan.(gg){ci}; qg=q{1}; qn=q{2};
      if size(qg,1)>=5, e=ellipseXY(qg,1); plot(ax,e(1,:),e(2,:),'-','Color',shGo{ci},'LineWidth',2.2); end
      if size(qn,1)>=5, e=ellipseXY(qn,1); plot(ax,e(1,:),e(2,:),'-','Color',shNo{ci},'LineWidth',2.2); end
    end
    xlim(ax,lg(1,:)); ylim(ax,lg(2,:)); set(ax,'XTickLabel',[],'YTickLabel',[]);
    if mi==1, title(ax,gg,'FontSize',11); end
    if gc==1, ylabel(ax,methods{mi,2},'FontWeight','bold','FontSize',11); end
  end
end
sgtitle(tl2,'D1-D6 ellipses overlaid (1-SD): Go=blue, NoGo=red, shade light\rightarrowdark = D1\rightarrowD6. Tight overlap = stereotyped across days');
exportgraphics(fig2,fullfile(outDir,'multiDR_ellipse_overlay.png'),'Resolution',150,'BackgroundColor','white');
try, exportgraphics(fig2,fullfile(outDir,'multiDR_ellipse_overlay.pdf'),'ContentType','image','BackgroundColor','white'); catch, end
try, savefig(fig2,fullfile(outDir,'multiDR_ellipse_overlay.fig')); catch, end
fprintf('Saved multiDR_ellipse_overlay.png to %s\n',outDir);

%% ---- covariance ellipse (k-SD) of 2D points ----
function xy=ellipseXY(P,k)
  mu=mean(P,1); C=cov(P); C=(C+C')/2; [V,D]=eig(C); D=max(D,0);
  t=linspace(0,2*pi,80); xy=(V*sqrt(D)*k)*[cos(t);sin(t)] + mu';
end

%% ---- two-pass alignment (same procedure as the Fig 2 cohort renders) ----
function A=twoPassAlignGroup(S,snapT)
  A=struct('cloudY',{},'cloudLab',{},'day',{},'goV',{},'noV',{});
  if isempty(S), return; end
  dim=size(S(1).Y,2);
  S2=struct('Cstack',{},'cloudY',{},'cloudLab',{},'day',{},'goV',{},'noV',{});
  for i=1:numel(S)
    twi=S(i).twi; twc=S(i).twc; Y=S(i).Y; lab=S(i).labels;
    uw=sort(unique(twi)); tC=arrayfun(@(w)mean(twc(twi==w),'omitnan'),uw); nW=numel(uw);
    tg=nan(nW,dim); tn=nan(nW,dim);
    for wi=1:nW
      g=Y(lab==S(i).goV & twi==uw(wi),:); n=Y(lab==S(i).noV & twi==uw(wi),:);
      if size(g,1)>=3, tg(wi,:)=mean(g,1); end
      if size(n,1)>=3, tn(wi,:)=mean(n,1); end
    end
    [~,wsi]=min(abs(tC-snapT)); sw=uw(wsi); m=twi==sw;
    S2(end+1)=struct('Cstack',[tg;tn],'cloudY',Y(m,:),'cloudLab',lab(m),'day',S(i).day,'goV',S(i).goV,'noV',S(i).noV); %#ok<AGROW>
  end
  days=unique([S2.day]); refByDay=containers.Map('KeyType','double','ValueType','any');
  for d=days
    di=find([S2.day]==d); if isempty(di), continue; end; nS=numel(di);
    if nS==1, refIdx=1; else
      DM=inf(nS);
      for a=1:nS, for b=1:nS, Ca=S2(di(a)).Cstack; Cb=S2(di(b)).Cstack; gd=all(isfinite(Ca),2)&all(isfinite(Cb),2); if sum(gd)>=6, DM(a,b)=mean(vecnorm(Ca(gd,:)-Cb(gd,:),2,2)); end, end, end
      [~,refIdx]=min(sum(DM,2,'omitnan'));
    end
    Cref=S2(di(refIdx)).Cstack;
    for a=1:nS
      Ci=S2(di(a)).Cstack; gd=all(isfinite(Cref),2)&all(isfinite(Ci),2);
      if sum(gd)>=6
        [~,~,tr]=procrustes(Cref(gd,:),Ci(gd,:),'Scaling',false,'Reflection','best');
        S2(di(a)).cloudY=tr.b*S2(di(a)).cloudY*tr.T+tr.c(1,:);
        S2(di(a)).Cstack=tr.b*Ci*tr.T+tr.c(1,:);
      end
    end
    refByDay(d)=S2(di(refIdx)).Cstack;
  end
  Cglob=refByDay(days(1));
  for d=days
    di=find([S2.day]==d); Crd=refByDay(d); gd=all(isfinite(Cglob),2)&all(isfinite(Crd),2);
    if sum(gd)>=6
      [~,~,trd]=procrustes(Cglob(gd,:),Crd(gd,:),'Scaling',false,'Reflection','best');
      for a=di, S2(a).cloudY=trd.b*S2(a).cloudY*trd.T+trd.c(1,:); end
    end
  end
  for i=1:numel(S2), A(end+1)=struct('cloudY',S2(i).cloudY,'cloudLab',S2(i).cloudLab,'day',S2(i).day,'goV',S2(i).goV,'noV',S2(i).noV); end %#ok<AGROW>
end
