% Export Figure 3 native-space geometry: clean stats file + long-format plot-data CSV.
% Reads geomTiming.mat (stats) + the session mat (per-session timecourses).
base   = fileparts(fileparts(mfilename('fullpath')));
sess   = load(fullfile(base,'cache','originalSpaceBasicGeometrySession.mat'));
gt     = load(fullfile(base,'stats','geomTiming.mat')); geomTiming = gt.geomTiming;
tvec   = sess.tvec(:)';
days   = sess.BehavCond;                 % {'D1'..'D6'}
groups = {'Naive','Trained'};
gng    = [90 45 30 25 20 15];

%% ---- clean stats file ----
fid = fopen(fullfile(base,'stats','native_geometry_stats.txt'),'w');
fprintf(fid,'FIGURE 3 - native-space geometry (centroid distance, radius)\n');
fprintf(fid,'Cluster-based permutation test (groupSeparationTest): %d perms, |Hedges g|>=%.1f, alpha=%.2f.\n', ...
    geomTiming.params.nPerms, geomTiming.params.minHedgesG, geomTiming.params.alpha);
fprintf(fid,'Naive vs Trained, per session. Onset = first sustained significant timepoint (s).\n');
fprintf(fid,'Tail: centDist right (trained>naive), radius left (trained<naive).\n');
fprintf(fid,'clusterP = fraction of %d permutations with max cluster mass >= observed (Maris & Oostenveld 2007), floored at 1e-4.\n\n', geomTiming.params.nPerms);
mm = {'centDist','Centroid distance'; 'radiusMix','Radius (mix)'};
for k=1:size(mm,1)
    mName=mm{k,1};
    fprintf(fid,'== %s ==\n', mm{k,2});
    fprintf(fid,'  Pooled (all days, n=%d naive / %d trained sessions): onset=%.3f s, clusterP=%.4g, nSigSust=%d/%d\n', ...
        geomTiming.(mName).nSessNaivePooled, geomTiming.(mName).nSessTrainPooled, ...
        geomTiming.(mName).onsetTimePooled, geomTiming.(mName).clusterPminPooled, ...
        sum(geomTiming.(mName).sigMaskSustPooled), numel(geomTiming.(mName).sigMaskSustPooled));
    fprintf(fid,'  Per day:  %-4s %-10s %-10s %-8s\n','day','onset(s)','clusterP','nSigSust');
    for d=1:numel(days)
        fprintf(fid,'            %-4s %-10s %-10.4g %-8d  (%d deg)\n', days{d}, ...
            num2str(geomTiming.(mName).onsetTime(d),'%.3f'), geomTiming.(mName).clusterPmin(d), ...
            geomTiming.(mName).nSigSust(d), gng(d));
    end
    fprintf(fid,'\n');
end
fclose(fid);
type(fullfile(base,'stats','native_geometry_stats.txt'));

%% ---- long-format plot-data CSV ----
metrics = {'centroidDistance', sess.centDist_bySess; ...
           'radius',           sess.radius_bySess; ...
           'radius_Go',        sess.radiusGo_bySess; ...
           'radius_NoGo',      sess.radiusNo_bySess; ...
           'centroidMag_Go',   sess.centMagGo_bySess; ...
           'centroidMag_NoGo', sess.centMagNo_bySess};
rows = {};
for mi=1:size(metrics,1)
    mName=metrics{mi,1}; C=metrics{mi,2};      % C is {2 x nDays}, each {1 x nSess} of [1 x nWin]
    for gi=1:2
        for d=1:numel(days)
            sc = C{gi,d}; if isempty(sc), continue; end
            for s=1:numel(sc)
                v=sc{s}; if isempty(v)||numel(v)~=numel(tvec), continue; end
                for w=1:numel(tvec)
                    rows(end+1,:)={mName, groups{gi}, days{d}, gng(d), s, tvec(w), v(w)}; %#ok<SAGROW>
                end
            end
        end
    end
end
T=cell2table(rows,'VariableNames',{'metric','group','day','gngAngle_deg','session','time_s','value'});
csvPath=fullfile(base,'data','native_geometry_plotdata_long.csv');
writetable(T,csvPath);
fprintf('\nWrote %d rows to %s\n', height(T), csvPath);
