%% multiDRembed.m  (Stage A: compute + cache embeddings for several DR methods)
% For each session, embed the SAME native data with: PaCMAP & UMATO (precomputed,
% read from *_PaCMAP_UMATO.mat), and classical MDS / Isomap / non-metric MDS
% (computed here from the forDR z-scored X). Points are stratified-subsampled per
% time window so trajectories + snapshot clouds are both covered. Resumable cache.

clear; clc;
% Baseline-z embeddings live in pacmap/cache/baselineZ (PaCMAP+UMATO from
% 02_make_pacmap_umato.py, with the forDR inputs co-located there).
codeDir=fileparts(mfilename('fullpath'));
baseRoot=getenv('MULTIDR_BASEROOT');
if isempty(baseRoot), baseRoot=fullfile(fileparts(fileparts(codeDir)),'pacmap','cache','baselineZ'); end
DBlist={'BatuhanADNaive','BatuhanAngDist'};
DBpretty=containers.Map(DBlist,{'Naive','Trained'});
days=1:6; goOri=45; nSub=1000; kIso=12; rng(0);
cacheFile=getenv('MULTIDR_CACHE');
if isempty(cacheFile), cacheFile=fullfile(fileparts(codeDir),'cache','multiDR_embeddings_baselinez.mat'); end
[cfDir,~,~]=fileparts(cacheFile); if ~exist(cfDir,'dir'), mkdir(cfDir); end

E=struct('group',{},'day',{},'animal',{},'goV',{},'noV',{},'labels',{},'twi',{},'twc',{}, ...
         'Ypac',{},'Yuma',{},'Ycmds',{},'Yiso',{},'Ynmds',{});
done=strings(0,1);
if exist(cacheFile,'file'), L=load(cacheFile); E=L.E; for s=1:numel(E), done(end+1)=E(s).group+"|"+E(s).animal+"|"+string(E(s).day); end; fprintf('Resuming: %d cached\n',numel(E)); end %#ok<SAGROW>

for dbi=1:numel(DBlist)
  DB=DBlist{dbi};
  for d=days
    pf=dir(fullfile(baseRoot,DB,sprintf('D%d',d),'*_PaCMAP_UMATO.mat'));
    for fi=1:numel(pf)
      base=erase(pf(fi).name,'_PaCMAP_UMATO.mat');
      key=string(DBpretty(DB))+"|"+string(strtok(base,'_'))+"_"+extractAfter(base,'_')+"|"+string(d); %#ok<NASGU>
      animal=strtok(erase(base,[DB '_']),'_');
      sk=string(DBpretty(DB))+"|"+string(animal)+"|"+string(d);
      if ismember(sk,done), continue; end
      P=load(fullfile(pf(fi).folder,pf(fi).name));
      forDR=fullfile(pf(fi).folder,[base '_TimeDR_forDR.mat']);
      if ~exist(forDR,'file'), continue; end
      F=load(forDR,'X','labels','time_window_idx','time_window_center');
      if ~isfield(P,'Y_pacmap')||~isfield(F,'X'), continue; end
      labels=double(P.labels(:)); twi=double(P.time_window_idx(:)); twc=double(P.time_window_center(:));
      X=double(F.X); if size(X,1)~=numel(labels), continue; end
      sv=unique(labels(isfinite(labels))); if numel(sv)~=2, continue; end
      if any(sv==goOri), goV=goOri; nv=sv(sv~=goOri); noV=nv(1); else, goV=sv(1); noV=sv(2); end
      % stratified subsample by window
      uw=unique(twi); per=max(1,round(nSub/numel(uw))); idx=[];
      for w=uw', r=find(twi==w); if numel(r)>per, r=r(randperm(numel(r),per)); end; idx=[idx;r]; end %#ok<AGROW>
      idx=sort(idx);
      lab=labels(idx); tw=twi(idx); tc=twc(idx); Xs=X(idx,:);
      % embeddings
      Ypac=double(P.Y_pacmap(idx,1:3));
      Yuma=double(P.Y_umato(idx,1:min(2,size(P.Y_umato,2))));
      Dmat=pdist(Xs);
      Ycmds=safeCmds(Dmat,3);
      Yiso=safeIso(Xs,kIso,3);
      Ynmds=safeNmds(Dmat,3);
      E(end+1)=struct('group',string(DBpretty(DB)),'day',d,'animal',string(animal),'goV',goV,'noV',noV, ...
        'labels',lab,'twi',tw,'twc',tc,'Ypac',Ypac,'Yuma',Yuma,'Ycmds',Ycmds,'Yiso',Yiso,'Ynmds',Ynmds); %#ok<SAGROW>
      done(end+1)=sk; %#ok<SAGROW>
      save(cacheFile,'E','-v7.3');
      fprintf('  %s D%d %s  [total %d]\n',char(DBpretty(DB)),d,animal,numel(E));
    end
  end
end
fprintf('Done. %d sessions cached at %s\n',numel(E),cacheFile);

function Y=safeCmds(Dmat,p)
  try, Y=cmdscale(Dmat,p); if size(Y,2)<p, Y=[Y zeros(size(Y,1),p-size(Y,2))]; end; catch, Y=[]; end
end
function Y=safeIso(Xs,k,p)
  try
    n=size(Xs,1); [idx,dd]=knnsearch(Xs,Xs,'K',k+1);
    ii=repmat((1:n)',1,k); jj=idx(:,2:end); ww=dd(:,2:end);
    A=sparse(ii(:),jj(:),ww(:),n,n); A=max(A,A'); G=graph(A);
    Dg=distances(G); if any(~isfinite(Dg(:))), Dg(~isfinite(Dg))=max(Dg(isfinite(Dg)))*1.5; end
    Y=cmdscale(Dg,p); if size(Y,2)<p, Y=[Y zeros(size(Y,1),p-size(Y,2))]; end
  catch, Y=[]; end
end
function Y=safeNmds(Dmat,p)
  try, Y=mdscale(Dmat,p,'Criterion','stress','Start','cmdscale','Options',statset('MaxIter',150,'TolFun',1e-3)); catch, Y=[]; end
end
