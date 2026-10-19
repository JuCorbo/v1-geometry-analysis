% Figure 3 dimensionality production driver.
%  - loads dimMetrics from the SVMPCAsession mat (authority: SVMPCATimeCourse_updated.m)
%  - renders the 1x4 summary figure (make_plotDimSVMSummaries.m)
%  - runs cluster-permutation group tests for ALL 3 metrics (D_PR, k80, tailSlope)
%  - writes stats file + long-format plot-data CSV
base    = fileparts(fileparts(mfilename('fullpath')));   % Figure3/dimensionality
addpath(fullfile(base,'code'));                          % make + local helpers
addpath(fullfile(fileparts(fileparts(base)),'_sharedWindowedData','code'));  % canonical groupSeparationTest
matPath = fullfile(base,'cache','dimMetrics_baseline012.mat');

fprintf('Loading dimMetrics ...\n');
S = load(matPath,'dimMetrics'); dimMetrics = S.dimMetrics; clear S;

%% ---- render the figure (uses dimMetrics in workspace; saves to panels\_raw) ----
run('make_plotDimSVMSummaries.m');

%% ---- build per-session Go+NoGo-averaged timecourse matrices per metric ----
dm_sz = ones(1,4); z=size(dimMetrics); dm_sz(1:numel(z))=z;
nDB=dm_sz(1); nDay=dm_sz(2); nSe=dm_sz(3); nRe=dm_sz(4);
eigSpecDay = 2;
% find tvec
tvec=[];
for dbi=1:nDB, for d=1:nDay, for s=1:nSe, for r=1:nRe
    dmc=dimMetrics(dbi,d,s,r);
    if ~isempty(dmc) && ~isempty(dmc.go), tvec=dmc.tvec(:)'; break; end
end, if ~isempty(tvec),break;end, end, if ~isempty(tvec),break;end, end, if ~isempty(tvec),break;end, end
nW=numel(tvec);
fields={'D_PR','k80','tailSlope'};

%% ---- stats: pooled (all days) + D2, for each metric ----
tails = struct('D_PR','left','k80','left','tailSlope','left'); % trained expected lower / steeper
fid=fopen(fullfile(base,'stats','dimensionality_stats.txt'),'w');
fprintf(fid,'FIGURE 3 - dimensionality (PCA eigenspectrum): D_PR (participation ratio), k80, tail slope\n');
fprintf(fid,'Compute: SVMPCATimeCourse_updated.m (dimMetrics). Go+NoGo averaged per window, then reps averaged within session (one value per session; n=42 naive / 60 trained).\n');
fprintf(fid,'Cluster-permutation group test (groupSeparationTest): 10000 perms, |Hedges g|>=0.2, alpha=0.05, tail=left (trained lower).\n');
fprintf(fid,'Onset = first sustained significant timepoint (s). dbi: 1=Naive, 2=Trained.\n');
fprintf(fid,'clusterP = fraction of 10,000 permutations with max cluster mass >= observed (Maris & Oostenveld 2007), floored at 1e-4.\n\n');
CSVrows={};
for fi=1:numel(fields)
    f=fields{fi};
    % matrices
    [Npool,~]=buildMat(dimMetrics,1,1:nDay,nSe,nRe,nW,{f}); [Tpool,~]=buildMat(dimMetrics,2,1:nDay,nSe,nRe,nW,{f});
    [Nd2,~]  =buildMat(dimMetrics,1,eigSpecDay,nSe,nRe,nW,{f}); [Td2,~]=buildMat(dimMetrics,2,eigSpecDay,nSe,nRe,nW,{f});
    sP=groupSeparationTest(Npool,Tpool,tvec,'tail',tails.(f),'nPerms',10000,'minSustained',1,'minHedgesG',0.2,'correction','cluster');
    sD=groupSeparationTest(Nd2, Td2, tvec,'tail',tails.(f),'nPerms',10000,'minSustained',1,'minHedgesG',0.2,'correction','cluster');
    % direction sanity at t~1s
    [~,w1]=min(abs(tvec-1.0)); dsign = mean(Tpool(:,w1),'omitnan')-mean(Npool(:,w1),'omitnan');
    fprintf(fid,'== %s ==\n',f);
    fprintf(fid,'  Pooled (n=%d naive / %d trained sessions): onset=%s s, clusterP=%.4g, nSigSust=%d/%d, (Trained-Naive)@1s=%.3f\n', ...
        size(Npool,1),size(Tpool,1), num2str(sP.onsetTime,'%.3f'), sP.clusterPmin, sum(sP.sigMaskSust), numel(sP.sigMaskSust), dsign);
    fprintf(fid,'  D2 only: onset=%s s, clusterP=%.4g, nSigSust=%d/%d\n\n', num2str(sD.onsetTime,'%.3f'), sD.clusterPmin, sum(sD.sigMaskSust), numel(sD.sigMaskSust));
    % CSV rows (per sess-rep, pooled over days label kept)
    for dbi=1:nDB
        for d=1:nDay
            [M,ids]=buildMat(dimMetrics,dbi,d,nSe,nRe,nW,{f});
            for i=1:size(M,1)
                for w=1:nW
                    CSVrows(end+1,:)={f, tern(dbi==1,'Naive','Trained'), sprintf('D%d',d), ids(i,2), tvec(w), M(i,w)}; %#ok<SAGROW>
                end
            end
        end
    end
end
fclose(fid); type(fullfile(base,'stats','dimensionality_stats.txt'));

T=cell2table(CSVrows,'VariableNames',{'metric','group','day','session','time_s','value'});
writetable(T,fullfile(base,'data','dimensionality_plotdata_long.csv'));
fprintf('\nWrote %d CSV rows.\n',height(T));

%% ---- local functions ----
function [M,ids]=buildMat(dimMetrics,dbi,days,nSe,nRe,nW,fields)
    % One row per SESSION: Go+NoGo averaged per window, then averaged across reps.
    f=fields{1}; M=[]; ids=[];
    for d=days(:)'
        for s=1:nSe
            acc=nan(nRe,nW);
            for r=1:nRe
                dmc=dimMetrics(dbi,d,s,r);
                if isempty(dmc)||isempty(dmc.go), continue; end
                nWloc=numel(dmc.go); row=nan(1,nW);
                for w=1:min(nWloc,nW)
                    row(w)=mean([dmc.go(w).(f), dmc.nogo(w).(f)],'omitnan');
                end
                acc(r,:)=row;
            end
            if all(isnan(acc(:))), continue; end
            M(end+1,:)=mean(acc,1,'omitnan'); ids(end+1,:)=[d s]; %#ok<AGROW>
        end
    end
end
function o=tern(c,a,b), if c, o=a; else, o=b; end, end
