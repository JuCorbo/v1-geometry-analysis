% PaCMAP-only overlaid-ellipses figure (Fig 2e: D1-D6, 1-SD Go/NoGo).
% baseRoot/output overridable from env (PAC_BASEROOT, PAC_OUT, PAC_TITLE) to run
% the baseline-z variant; defaults render the main figure.
codeDir=fileparts(mfilename('fullpath'));
baseRoot=getenv('PAC_BASEROOT'); outPng=getenv('PAC_OUT'); ttl=getenv('PAC_TITLE');
if isempty(baseRoot), baseRoot=fullfile(codeDir,'..','cache','main_pacmap_bwgz_pm080','baseline_seed0'); end
if isempty(outPng),   outPng=fullfile(codeDir,'..','panels','main','e_overlay_D1toD6_ellipses.png'); end
if isempty(ttl),      ttl='Main (global-z + baseline subtraction)'; end
DBlist={'BatuhanADNaive','BatuhanAngDist'}; DBpretty=containers.Map(DBlist,{'Naive','Trained'});
days=1:6; goOri=45; nSub=1000; snapT=2.5; rng(0); daysToPlot=1:6; nD=numel(daysToPlot);
colGo=[0.20 0.40 0.80]; colNo=[0.80 0.20 0.20];

% ---- build E (PaCMAP only) ----
E=struct('group',{},'day',{},'labels',{},'twi',{},'twc',{},'goV',{},'noV',{},'Ypac',{});
for dbi=1:numel(DBlist)
  DB=DBlist{dbi};
  for d=days
    pf=dir(fullfile(baseRoot,DB,sprintf('D%d',d),'*_PaCMAP_UMATO.mat'));
    for fi=1:numel(pf)
      P=load(fullfile(pf(fi).folder,pf(fi).name));
      if ~isfield(P,'Y_pacmap'), continue; end
      labels=double(P.labels(:)); twi=double(P.time_window_idx(:)); twc=double(P.time_window_center(:));
      sv=unique(labels(isfinite(labels))); if numel(sv)~=2, continue; end
      if any(sv==goOri), goV=goOri; nv=sv(sv~=goOri); noV=nv(1); else, goV=sv(1); noV=sv(2); end
      uw=unique(twi); per=max(1,round(nSub/numel(uw))); idx=[];
      for w=uw', r=find(twi==w); if numel(r)>per, r=r(randperm(numel(r),per)); end; idx=[idx;r]; end %#ok<AGROW>
      idx=sort(idx);
      E(end+1)=struct('group',string(DBpretty(DB)),'day',d,'labels',labels(idx),'twi',twi(idx), ...
        'twc',twc(idx),'goV',goV,'noV',noV,'Ypac',double(P.Y_pacmap(idx,1:3))); %#ok<AGROW>
    end
  end
end
grp=string({E.group})';

% ---- Phase 1: align + discriminant projection (PaCMAP) ----
GA=struct(); PgAll=[]; PnAll=[];
for gg=["Naive","Trained"]
  S=struct('Y',{},'labels',{},'twi',{},'twc',{},'goV',{},'noV',{},'day',{});
  for j=find(grp==gg)'
    S(end+1)=struct('Y',E(j).Ypac,'labels',E(j).labels,'twi',E(j).twi,'twc',E(j).twc,'goV',E(j).goV,'noV',E(j).noV,'day',E(j).day); %#ok<AGROW>
  end
  A=twoPassAlignGroup(S,snapT); GA.(gg)=A;
  for a=1:numel(A)
    PgAll=[PgAll; A(a).cloudY(A(a).cloudLab==A(a).goV,:)]; %#ok<AGROW>
    PnAll=[PnAll; A(a).cloudY(A(a).cloudLab==A(a).noV,:)]; %#ok<AGROW>
  end
end
u1=(mean(PnAll,1)-mean(PgAll,1))'; u1=u1/norm(u1);
allP=[PgAll;PnAll]; res=allP-(allP*u1)*u1'; cf=pca(res); u2=cf(:,1);
W=[u1 u2]; proj=@(P)P*W;
pan=struct(); lim=struct();
for gg=["Naive","Trained"]
  A=GA.(gg); gAcc=[];
  for ci=1:nD
    d=daysToPlot(ci); Pg=[]; Pn=[];
    for a=1:numel(A)
      if A(a).day~=d, continue; end
      Pg=[Pg; A(a).cloudY(A(a).cloudLab==A(a).goV,:)]; Pn=[Pn; A(a).cloudY(A(a).cloudLab==A(a).noV,:)]; %#ok<AGROW>
    end
    qg=[]; qn=[]; if ~isempty(Pg), qg=proj(Pg); end; if ~isempty(Pn), qn=proj(Pn); end
    pan.(gg){ci}={qg,qn}; gAcc=[gAcc; qg; qn]; %#ok<AGROW>
  end
  xl=prctile(gAcc(:,1),[1.5 98.5]); yl=prctile(gAcc(:,2),[1.5 98.5]);
  xl=xl+0.12*[-1 1]*diff(xl); yl=yl+0.12*[-1 1]*diff(yl); lim.(gg)=[xl; yl];
end

% ---- render overlay (Naive | Trained) ----
blend=@(c,f)(1-f)*c+f*[1 1 1];
shGo=cell(1,nD); shNo=cell(1,nD);
for d=1:nD, fr=0.62*(nD-d)/max(nD-1,1); shGo{d}=blend(colGo,fr); shNo{d}=blend(colNo,fr); end
fig=figure('Color','w','Position',[60 60 900 460]);
tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for gg=["Naive","Trained"]
  ax=nexttile(tl); hold(ax,'on'); box(ax,'on'); grid(ax,'on');
  for ci=1:nD
    q=pan.(gg){ci}; qg=q{1}; qn=q{2};
    if size(qg,1)>=5, e=ellipseXY(qg,1); plot(ax,e(1,:),e(2,:),'-','Color',shGo{ci},'LineWidth',2.2); end
    if size(qn,1)>=5, e=ellipseXY(qn,1); plot(ax,e(1,:),e(2,:),'-','Color',shNo{ci},'LineWidth',2.2); end
  end
  xlim(ax,lim.(gg)(1,:)); ylim(ax,lim.(gg)(2,:)); set(ax,'XTickLabel',[],'YTickLabel',[]);
  title(ax,gg,'FontSize',12);
end
sgtitle(tl,sprintf('%s | PaCMAP D1-D6 ellipses (Go blue, NoGo red; light->dark = D1->D6)',ttl));
exportgraphics(fig,outPng,'Resolution',150,'BackgroundColor','white');
exportgraphics(fig,strrep(outPng,'.png','.pdf'),'ContentType','vector','BackgroundColor','white');
fprintf('saved %s\n',outPng);

function xy=ellipseXY(P,k)
  mu=mean(P,1); C=cov(P); C=(C+C')/2; [V,D]=eig(C); D=max(D,0);
  t=linspace(0,2*pi,80); xy=(V*sqrt(D)*k)*[cos(t);sin(t)] + mu';
end
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
