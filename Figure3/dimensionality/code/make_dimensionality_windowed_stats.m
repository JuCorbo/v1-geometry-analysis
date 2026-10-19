%% make_dimensionality_windowed_stats.m
% Windowed group statistics for the dimensionality measures, in the SAME
% convention as the other measures (see nativeSpaceGeometry/stats/peak*_byDay):
%   - scalar windows (baseline t<0, 1-3 s, 0-3 s): per-session means (Go/NoGo and
%     rounds averaged), naive vs trained Wilcoxon rank-sum + Hedges g (pooled SD).
%   - baseline -> stimulus within group: paired Wilcoxon signed-rank + Hedges
%     g_av (RMS of the two within-condition SDs), per Statistics section.
%   - pooled-classes "overall" PR (cache/pooledSpectralPR.mat, from
%     pull_pooledSpectralPR.m): cluster-permutation timecourse test + the same
%     scalar windows.
% Writes stats/dimensionality_windowed_stats.txt.
clear; clc;
codeDir = fileparts(mfilename('fullpath'));
aDir    = fileparts(codeDir); figDir = fileparts(fileparts(aDir));
addpath(fullfile(figDir,'_sharedWindowedData','code'));   % groupSeparationTest
L = load(fullfile(aDir,'cache','dimMetrics_baseline012.mat'));
dimMetrics = L.dimMetrics;
P = load(fullfile(aDir,'cache','pooledSpectralPR.mat'));

dm_sz = ones(1,4); z = size(dimMetrics); dm_sz(1:numel(z)) = z;
tvec = [];
for dbi=1:dm_sz(1), for d=1:dm_sz(2), for s=1:dm_sz(3), for r=1:dm_sz(4)
    dmc = dimMetrics(dbi,d,s,r);
    if ~isempty(dmc) && ~isempty(dmc.go), tvec = dmc.tvec(:)'; break; end
end, if ~isempty(tvec),break;end, end, if ~isempty(tvec),break;end, end, if ~isempty(tvec),break;end, end
nW = numel(tvec);
wins = struct('name',{'base','1-3s','0-3s'}, 'mask',{tvec<0, tvec>=1&tvec<=3, tvec>=0&tvec<=3});

fields = {'D_PR','k80','tailSlope'};
fid = fopen(fullfile(aDir,'stats','dimensionality_windowed_stats.txt'),'w');
fprintf(fid,'FIGURE 3 - dimensionality: windowed group statistics\n');
fprintf(fid,'One value per session: Go+NoGo averaged per window, rounds averaged, mean over window; pooled across days.\n');
fprintf(fid,'Between-group: Wilcoxon rank-sum, Hedges g (pooled SD, small-sample corrected), mean +/- SEM reported.\n');
fprintf(fid,'Within-group baseline->stimulus: Wilcoxon signed-rank, Hedges g_av (RMS of the two SDs).\n');
fprintf(fid,'Timecourse tests (onset, cluster p) live in dimensionality_stats.txt.\n\n');

%% ---- within-class measures (dimMetrics) ----
for fi = 1:numel(fields)
    f = fields{fi};
    W = cell(2, numel(wins));                       % per-session window means
    for dbi = 1:2
        TC = sessionTCs_dm(dimMetrics, dbi, dm_sz, nW, f);
        for wi = 1:numel(wins)
            W{dbi,wi} = mean(TC(:, wins(wi).mask), 2, 'omitnan');
        end
    end
    fprintf(fid,'== %s ==\n', f);
    for wi = 2:3                                    % 1-3s, 0-3s between-group
        writeUnpaired(fid, wins(wi).name, W{1,wi}, W{2,wi});
    end
    writeUnpaired(fid, 'baseline (t<0)', W{1,1}, W{2,1});
    for dbi = 1:2
        writePaired(fid, sprintf('%s baseline->1-3s', tern(dbi==1,'Naive  ','Trained')), W{dbi,1}, W{dbi,2});
    end
    fprintf(fid,'\n');
end

%% ---- pooled-classes spectral PR ----
fprintf(fid,'== Pooled-classes ("overall") D_PR ==  (pull_pooledSpectralPR.m: PCA on Go+NoGo pooled,\n');
fprintf(fid,'   pooled-mean centred; the spectrum contains both within-class variability and between-class separation)\n');
TCp = cell(2,1);
for dbi = 1:2
    M = [];
    for d = 1:size(P.PRpool,2)
        for s = 1:numel(P.PRpool{dbi,d})
            M(end+1,:) = mean(P.PRpool{dbi,d}{s}, 1, 'omitnan'); %#ok<SAGROW>
        end
    end
    TCp{dbi} = M;
end
sPool = groupSeparationTest(TCp{1}, TCp{2}, tvec, 'tail','left', 'nPerms',10000, ...
    'minSustained',1, 'minHedgesG',0.2, 'correction','cluster');
fprintf(fid,'  Timecourse (cluster permutation, 10000 perms, tail left): onset=%s s, clusterP=%.4g, nSigSust=%d/%d\n', ...
    num2str(sPool.onsetTime,'%.3f'), sPool.clusterPmin, sum(sPool.sigMaskSust), numel(sPool.sigMaskSust));
for wi = 1:numel(wins)
    vN = mean(TCp{1}(:, wins(wi).mask), 2, 'omitnan');
    vT = mean(TCp{2}(:, wins(wi).mask), 2, 'omitnan');
    writeUnpaired(fid, wins(wi).name, vN, vT);
end
fclose(fid);
type(fullfile(aDir,'stats','dimensionality_windowed_stats.txt'));

%% ---- helpers ----
function TC = sessionTCs_dm(dimMetrics, dbi, dm_sz, nW, f)
    TC = [];
    for d = 1:dm_sz(2)
        for s = 1:dm_sz(3)
            acc = nan(dm_sz(4), nW);
            for r = 1:dm_sz(4)
                dmc = dimMetrics(dbi,d,s,r);
                if isempty(dmc) || isempty(dmc.go), continue; end
                row = nan(1,nW);
                for w = 1:min(numel(dmc.go), nW)
                    row(w) = mean([dmc.go(w).(f), dmc.nogo(w).(f)], 'omitnan');
                end
                acc(r,:) = row;
            end
            if all(isnan(acc(:))), continue; end
            TC(end+1,:) = mean(acc, 1, 'omitnan'); %#ok<AGROW>
        end
    end
end
function writeUnpaired(fid, lab, vN, vT)
    vN = vN(isfinite(vN)); vT = vT(isfinite(vT));
    n1 = numel(vN); n2 = numel(vT);
    sp = sqrt(((n1-1)*var(vN) + (n2-1)*var(vT)) / (n1+n2-2));
    g  = (1 - 3/(4*(n1+n2)-9)) * (mean(vT)-mean(vN)) / sp;
    fprintf(fid,'  %-16s: Naive %.4g +/- %.2g (n=%d) | Trained %.4g +/- %.2g (n=%d) | rank-sum p=%.4g | Hedges g=%+.2f\n', ...
        lab, mean(vN), std(vN)/sqrt(n1), n1, mean(vT), std(vT)/sqrt(n2), n2, ranksum(vN,vT), g);
end
function writePaired(fid, lab, vA, vB)
    ok = isfinite(vA) & isfinite(vB); vA = vA(ok); vB = vB(ok); n = numel(vA);
    sav = sqrt((var(vA)+var(vB))/2);
    gav = (1 - 3/(4*(n-1)-1)) * (mean(vB)-mean(vA)) / sav;
    fprintf(fid,'  %-24s: %.4g -> %.4g (n=%d) | signed-rank p=%.4g | Hedges g_av=%+.2f\n', ...
        lab, mean(vA), mean(vB), n, signrank(vA,vB), gav);
end
function o = tern(c,a,b), if c, o=a; else, o=b; end, end
