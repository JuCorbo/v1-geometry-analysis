%% writeTimingStats.m — dump the reported timing statistics (Figure 6 a/b) to stats/
%  Reads ../data/onsetCI.mat (neural divergence onsets + bootstrap CIs) and
%  ../data/DecisionUB_ALL_BatuhanAngDist.csv (behavioural decision proxies).
clear; clc;
codeDir=fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
anaDir=fileparts(codeDir);
R=load(fullfile(anaDir,'data','onsetCI.mat')); res=R.results;
T=readtable(fullfile(anaDir,'data','DecisionUB_ALL_BatuhanAngDist.csv'));
outFile=fullfile(anaDir,'stats','timing_stats.txt');
fid=fopen(outFile,'w');

fprintf(fid,'FIGURE 6 timing thread — reported statistics\n');
fprintf(fid,'Convention: neural divergence onsets reported at the window CENTRE (tvec_centre), 0.12 preprocessing.\n');
fprintf(fid,'Onsets = pooled naive-vs-trained cluster-permutation onset; CI = 95%% over 1000 session bootstraps.\n\n');

fprintf(fid,'NEURAL DIVERGENCE ONSETS (s from stimulus onset)\n');
measOrder={'dimDPR','radiusMix','capacity','centDist'};
lbl=struct('dimDPR','Dimensionality (D_PR)','radiusMix','Manifold radius','capacity','Manifold capacity','centDist','Centroid distance');
for i=1:numel(measOrder)
  m=measOrder{i}; if ~isfield(res,m), continue; end; s=res.(m);
  fprintf(fid,'  %-22s onset=%.3f  95%% CI [%.3f, %.3f]  bootMed=%.3f  detected=%.0f%%  (nN=%d nT=%d, tail=%s)\n', ...
    lbl.(m), s.onsetPoint, s.ci(1), s.ci(2), s.median, 100*s.detRate, s.nNaive, s.nTrained, s.tail);
end

fprintf(fid,'\nBEHAVIOURAL DECISION PROXIES (trained; s from stimulus onset)\n');
vn=T.Properties.VariableNames;
dayCol='BehavCond'; if ~ismember(dayCol,vn), dayCol=vn{find(contains(lower(vn),'cond'),1)}; end
days=unique(string(T.(dayCol)));
proxCols={'tUpperBound','medCR_last'};
proxName={'Go/NoGo lick divergence','Last correct-rejection lick'};
for p=1:numel(proxCols)
  c=proxCols{p}; if ~ismember(c,vn), continue; end
  v=double(T.(c)); g=string(T.(dayCol));
  ok=isfinite(v);
  grp=findgroups(g(ok)); pKW=kruskalwallis(v(ok),grp,'off');
  fprintf(fid,'  %-30s median=%.3f  Kruskal-Wallis day effect p=%.3g\n', proxName{p}, median(v(ok)), pKW);
  for d=1:numel(days)
    sel=ok & g==days(d);
    if any(sel), fprintf(fid,'      %s median=%.3f (n=%d)\n', days(d), median(v(sel)), sum(sel)); end
  end
end
fclose(fid);
fprintf('wrote %s\n', outFile);
type(outFile);
