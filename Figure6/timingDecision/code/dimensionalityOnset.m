%% Dimensionality (participation ratio) naive-vs-trained divergence onset
% Extract D_PR timecourse per session from dimMetrics, run groupSeparationTest
% with the SAME settings used for radius/centDist (tail='left' since training
% is expected to REDUCE dimensionality). Add onset to the timing comparison.
clear; clc;

codeDir = fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
anaDir = fileparts(codeDir);
outDir = codeDir;
figRoot = fileparts(fileparts(anaDir));
addpath(fullfile(figRoot,'_sharedWindowedData','code'));   % groupSeparationTest
srcMat = fullfile(figRoot,'Figure3','dimensionality','cache','dimMetrics_baseline012.mat');  % 0.12 uniformized dim

S = load(srcMat, 'dimMetrics', 'DBlist', 'BehavCond', 'tvec');
dimMetrics = S.dimMetrics; tvec = S.tvec(:)';
days = S.BehavCond; nDay = numel(days);
fprintf('DBlist = %s\n', strjoin(cellstr(string(S.DBlist)), ' , '));

[nDB, nExp, nSess, nRep] = size(dimMetrics);

% identify naive vs trained by DB name
dbNames = cellstr(string(S.DBlist));
naiveDB = find(contains(lower(dbNames),'naive'), 1);
trainDB = find(~contains(lower(dbNames),'naive'), 1);
fprintf('naive DB idx=%d (%s), trained DB idx=%d (%s)\n', ...
    naiveDB, dbNames{naiveDB}, trainDB, dbNames{trainDB});

metric = 'D_PR';   % participation ratio = effective dimensionality
onsetOffsetSec = 1*0.0646;  % 1*dt: report at window centre (tvec_centre)

% ---- collect per-session mean timecourse (mean over reps, mean of Go & NoGo) ----
getTC = @(sc) arrayfun(@(x) x.(metric), sc(:))';   % 1 x nWin from a 26x1 struct
collect = @(dbi,ex) local_collect(dimMetrics, dbi, ex, nSess, nRep, getTC, numel(tvec));

dimN = cell(1,nDay); dimT = cell(1,nDay);
for ex = 1:nDay
    dimN{ex} = collect(naiveDB, ex);   % [nSessNaive x nWin]
    dimT{ex} = collect(trainDB, ex);
end

%% ---- per-day + pooled onset via groupSeparationTest ----
dimTiming.measure='D_PR'; dimTiming.tvec=tvec; dimTiming.days=days;
dimTiming.onsetTime=nan(1,nDay); dimTiming.nSigSust=zeros(1,nDay);
allN=[]; allT=[];
fprintf('\nPER-DAY dimensionality (D_PR) divergence (tail=left: trained<naive):\n');
fprintf('%-4s | nN nT | onset(s) | nSigSust | meanG(post)\n','Day');
for ex=1:nDay
    dN=dimN{ex}; dT=dimT{ex};
    if isempty(dN)||isempty(dT), fprintf('%-4s | skip\n',days{ex}); continue; end
    r = groupSeparationTest(dN, dT, tvec, 'tail','left','correction','cluster', ...
        'nPerms',10000,'alpha',0.05,'minSustained',1,'minHedgesG',0.2,'tStart',0, ...
        'onsetOffsetSec',onsetOffsetSec,'measureName',sprintf('D_PR_%s',days{ex}),'verbose',false);
    dimTiming.onsetTime(ex)=r.onsetTime; dimTiming.nSigSust(ex)=sum(r.sigMaskSust);
    dimTiming.hedgesG{ex}=r.hedgesG;
    fprintf('%-4s | %2d %2d | %8s | %8d | %.2f\n', days{ex}, size(dN,1), size(dT,1), ...
        num2str(r.onsetTime,'%.2f'), sum(r.sigMaskSust), mean(r.hedgesG(isfinite(r.hedgesG))));
    allN=[allN; dN]; allT=[allT; dT]; %#ok<AGROW>
end
rp = groupSeparationTest(allN, allT, tvec, 'tail','left','correction','cluster', ...
    'nPerms',10000,'alpha',0.05,'minSustained',1,'minHedgesG',0.2,'tStart',0, ...
    'onsetOffsetSec',onsetOffsetSec,'measureName','D_PR_pooled','verbose',false);
dimTiming.onsetPooled=rp.onsetTime; dimTiming.statsPooled=rp;
fprintf('POOLED (%d naive, %d trained): onset = %.3f s, nSigSust=%d, |maxG|=%.2f\n', ...
    size(allN,1), size(allT,1), rp.onsetTime, sum(rp.sigMaskSust), max(abs(rp.hedgesG(isfinite(rp.hedgesG)))));

% mean D_PR level (context: does training reduce dimensionality?)
mnN = mean(allN(:,tvec>=0.5 & tvec<=2.5),'all','omitnan');
mnT = mean(allT(:,tvec>=0.5 & tvec<=2.5),'all','omitnan');
fprintf('mean D_PR in [0.5,2.5]s: naive=%.2f  trained=%.2f  (trained/naive=%.2f)\n', mnN, mnT, mnT/mnN);

save(fullfile(anaDir,'cache','dimTiming.mat'),'dimTiming','dimN','dimT','tvec','days');
fprintf('\nSaved dimTiming.mat in %s\n', outDir);

function M = local_collect(dimMetrics, dbi, ex, nSess, nRep, getTC, nWin)
    M = [];
    for s=1:nSess
        acc=[];
        for r=1:nRep
            el = dimMetrics(dbi,ex,s,r);
            if isempty(el) || ~isfield(el,'go') || isempty(el.go), continue; end
            try
                g = getTC(el.go); n = getTC(el.nogo);
            catch, continue; end
            if numel(g)~=nWin, continue; end
            acc(end+1,:) = mean([g;n],1,'omitnan'); %#ok<AGROW>  % mean Go & NoGo
        end
        if ~isempty(acc)
            M(end+1,:) = mean(acc,1,'omitnan'); %#ok<AGROW>       % mean over reps
        end
    end
end
