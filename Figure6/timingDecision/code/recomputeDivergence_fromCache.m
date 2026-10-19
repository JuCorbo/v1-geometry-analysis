%% recomputeDivergence_fromCache.m
%  Builder of the Go/NoGo lick-divergence upper bound (tUpperBound) from the
%  cached full lick traces (lickTracesFull/). Compares against the shipped
%  DecisionUB CSV, then writes a *_fromCache.csv (does NOT clobber the original).
clear; clc;
here    = fileparts(mfilename('fullpath')); if isempty(here), here=pwd; end
anaDir  = fileparts(here);
csvPath = fullfile(anaDir,'data','DecisionUB_ALL_BatuhanAngDist.csv');
lickCacheRoot=fullfile(fileparts(fileparts(anaDir)),'Figure5','behavioralControls','cache');
lickDirs= {fullfile(lickCacheRoot,'lickTracesFull'), fullfile(lickCacheRoot,'lickTracesFull_D2')};
db='BatuhanAngDist'; goVal=45;
noGoByDay=containers.Map({'D1','D2','D3','D4','D5','D6'},{135,90,75,70,65,60});
stimOnsetSec=1.5; plotWinSec=[-1 6]; doSmooth=true; smoothWinSec=0.15;
threshMethod='kmeans2'; maxThreshPts=200000; clipPrct=[0.5 99.5]; %#ok<NASGU>
nBoot=300; ciAlpha=0.05; deltaMin=0.10; holdSec=0.20; ubSearch=[0 6];
rng(1);

% map Experiment -> cached lick-traces file (across both cache folders)
fileMap=containers.Map;
for d=1:numel(lickDirs)
    if ~isfolder(lickDirs{d}), continue; end
    ff=dir(fullfile(lickDirs{d}, sprintf('*_%s_lickTraces.mat',db)));
    for f=1:numel(ff)
        id=extractBefore(ff(f).name, sprintf('_%s_lickTraces.mat',db));
        if ~isempty(id) && ~isKey(fileMap,id), fileMap(id)=fullfile(ff(f).folder,ff(f).name); end
    end
end
fprintf('cached lick-trace sessions available: %d\n', fileMap.Count);

T=readtable(csvPath);
newUB=nan(height(T),1); nMiss=0;
for r=1:height(T)
    day=char(string(T.BehavCond(r))); expID=char(string(T.Experiment(r))); noGoVal=noGoByDay(day);
    if ~isKey(fileMap,expID), nMiss=nMiss+1; continue; end
    L=load(fileMap(expID),'lickTrace','Visual_Stim','Lick_dt','Trial');
    [~,ord]=sort(double(L.Trial)); vs=double(L.Visual_Stim(ord));      % match SQL ORDER BY Trial
    traces=L.lickTrace(ord);                                          % already-decoded cell
    keep = vs==goVal | vs==noGoVal; vs=vs(keep); traces=traces(keep);
    lickDt=median(double(L.Lick_dt(isfinite(double(L.Lick_dt)))));
    thr=computeSessionLickThreshold(traces,lickDt,plotWinSec,stimOnsetSec,doSmooth,smoothWinSec,threshMethod,maxThreshPts,clipPrct);
    if ~isfinite(thr), thr=computeFallbackMedianThreshold(traces,lickDt,plotWinSec,stimOnsetSec,doSmooth,smoothWinSec); end
    isGo=vs==goVal; isNoGo=vs==noGoVal;
    [tVec,goOn]=tracesToLickOn(traces(isGo),lickDt,stimOnsetSec,plotWinSec,doSmooth,smoothWinSec,thr);
    [~,noOn]  =tracesToLickOn(traces(isNoGo),lickDt,stimOnsetSec,plotWinSec,doSmooth,smoothWinSec,thr);
    if isempty(goOn)||isempty(noOn)||size(goOn,1)<5||size(noOn,1)<5, continue; end
    [~,~,delta,deltaLo]=bootstrapDeltaCI(goOn,noOn,nBoot,ciAlpha);
    newUB(r)=findFirstSustained(tVec,delta,deltaLo,deltaMin,holdSec,ubSearch);
end

oldUB=T.tUpperBound;
d=newUB-oldUB;
fprintf('sessions missing a cache file: %d / %d\n', nMiss, height(T));
fprintf('cache vs SQL tUpperBound: median|Δ|=%.4f s, max|Δ|=%.4f s, agree(<50ms)=%d/%d\n', ...
    median(abs(d),'omitnan'), max(abs(d)), sum(abs(d)<0.05), sum(isfinite(d)));

outCsv=strrep(csvPath,'.csv','_fromCache.csv');
Tout=T; Tout.tUpperBound=newUB; writetable(Tout,outCsv);
fprintf('SQL-free tUpperBound written to %s\n', outCsv);

%% helpers (identical to recomputeDivergence.m)
function [tVec,lickOn]=tracesToLickOn(traces,dt,stimOnsetSec,win,doSmooth,smoothWinSec,thr)
    tVec=[];lickOn=[]; if isempty(traces)||~(isfinite(dt)&&dt>0), return; end
    lens=cellfun(@numel,traces); Lmax=max(lens); if Lmax<2, return; end
    M=NaN(numel(traces),Lmax);
    for k=1:numel(traces), v=traces{k}; if isempty(v), continue; end; v=double(v(:)');
        if doSmooth&&smoothWinSec>0, v=movmean(v,max(1,round(smoothWinSec/dt)),'Endpoints','shrink'); end
        M(k,1:numel(v))=v; end
    tAll=(0:Lmax-1)*dt - stimOnsetSec; idx=tAll>=win(1)&tAll<=win(2); if ~any(idx), return; end
    tVec=tAll(idx); lickOn=M(:,idx)>thr;
end
function [pGo,pNo,delta,deltaLo]=bootstrapDeltaCI(goOn,noOn,nBoot,alpha)
    tLen=size(goOn,2); pGo=mean(goOn,1,'omitnan'); pNo=mean(noOn,1,'omitnan'); delta=pGo-pNo;
    bd=NaN(nBoot,tLen); nG=size(goOn,1); nN=size(noOn,1);
    for b=1:nBoot, bd(b,:)=mean(goOn(randi(nG,nG,1),:),1,'omitnan')-mean(noOn(randi(nN,nN,1),:),1,'omitnan'); end
    deltaLo=prctile(bd,100*(alpha/2),1);
end
function tStar=findFirstSustained(tVec,delta,deltaLo,deltaMin,holdSec,searchWinSec)
    tStar=NaN; if isempty(tVec), return; end
    idxSearch=tVec>=searchWinSec(1)&tVec<=searchWinSec(2); if ~any(idxSearch), return; end
    ok=deltaLo>0; if isfinite(deltaMin)&&deltaMin>0, ok=ok&(delta>deltaMin); end
    ok=ok&idxSearch; dt=median(diff(tVec)); nHold=max(1,round(holdSec/dt));
    okv=double(ok(:))'; dd=diff([0 okv 0]); s=find(dd==1); e=find(dd==-1)-1; g=find((e-s+1)>=nHold,1,'first');
    if ~isempty(g), tStar=tVec(s(g)); end
end
function thr=computeSessionLickThreshold(traces,dt,win,stimOnsetSec,doSmooth,smoothWinSec,method,maxPts,clipPrct) %#ok<INUSD>
    thr=NaN; vals=[];
    for k=1:numel(traces), v=traces{k}; if isempty(v), continue; end; v=double(v(:)');
        if doSmooth&&smoothWinSec>0, v=movmean(v,max(1,round(smoothWinSec/dt)),'Endpoints','shrink'); end
        t=(0:numel(v)-1)*dt - stimOnsetSec; idx=t>=win(1)&t<=win(2); if any(idx), vals=[vals v(idx)]; end; end
    vals=vals(isfinite(vals)); if numel(vals)<1000, return; end
    if numel(vals)>maxPts, vals=vals(randperm(numel(vals),maxPts)); end
    lo=prctile(vals,clipPrct(1)); hi=prctile(vals,clipPrct(2)); vals=vals(vals>=lo&vals<=hi); if numel(vals)<1000, return; end
    try, [~,C]=kmeans(vals(:),2,'Replicates',5,'MaxIter',200,'Display','off'); C=sort(C(:)); if all(isfinite(C))&&abs(C(2)-C(1))>eps, thr=mean(C); end; catch, end
end
function thr=computeFallbackMedianThreshold(traces,dt,win,stimOnsetSec,doSmooth,smoothWinSec)
    vals=[]; for k=1:numel(traces), v=traces{k}; if isempty(v), continue; end; v=double(v(:)');
        if doSmooth&&smoothWinSec>0, v=movmean(v,max(1,round(smoothWinSec/dt)),'Endpoints','shrink'); end
        t=(0:numel(v)-1)*dt - stimOnsetSec; idx=t>=win(1)&t<=win(2); if any(idx), vals=[vals v(idx)]; end; end
    thr=median(vals(isfinite(vals)),'omitnan');
end
