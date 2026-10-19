% Figure 2d theta stats + per-day summary CSV, from the main (global-z+sub, 0.8.0, seed0) run.
% Reads theta_bwgz_s0.csv (exported by s03); writes the deliverable CSVs to
% pacmap/data and theta_stats.txt to pacmap/stats. Run s03 first.
SP=fullfile(fileparts(mfilename('fullpath')),'..','data');   % where s03 exports theta_bwgz_s0.csv (or set F2_THETA_OUT)
codeDir=fileparts(mfilename('fullpath')); aDir=fileparts(codeDir);
dataDir=fullfile(aDir,'data'); statDir=fullfile(aDir,'stats');
if ~exist(dataDir,'dir'), mkdir(dataDir); end
if ~exist(statDir,'dir'), mkdir(statDir); end

T=readtable(fullfile(SP,'theta_bwgz_s0.csv'),'TextType','string');   % sessKey,day,group,theta
day=double(T.day); grp=string(T.group); th=double(T.theta);
gngMap=[90 45 30 25 20 15]; gng=gngMap(day)';

% per-session CSV
writetable(table(day,gng,grp,th,'VariableNames',{'day','gngAngle_deg','group','trajAngle_theta_deg'}), fullfile(dataDir,'theta_perSession.csv'));
% per-day x group summary
rows={};
for d=1:6, for g=["Naive","Trained"]
  v=th(day==d & grp==g); v=v(~isnan(v)); q=prctile(v,[25 50 75]);
  rows(end+1,:)={d,gngMap(d),char(g),numel(v),q(2),q(1),q(3),mean(v),std(v)}; %#ok<SAGROW>
end, end
writetable(cell2table(rows,'VariableNames',{'day','gngAngle_deg','group','n','median','q25','q75','mean','sd'}), fullfile(dataDir,'theta_perDay_groupSummary.csv'));

% stats
p_kw_group=kruskalwallis(th,grp,'off'); p_kw_day=kruskalwallis(th,day,'off');
p_kw_day_N=kruskalwallis(th(grp=="Naive"),day(grp=="Naive"),'off');
p_kw_day_T=kruskalwallis(th(grp=="Trained"),day(grp=="Trained"),'off');
fid=fopen(fullfile(statDir,'theta_stats.txt'),'w');
fprintf(fid,'FIGURE 2 - PaCMAP Go/NoGo trajectory angle (theta)\nGenerated %s\n\n',datestr(now));
fprintf(fid,'Input: main PaCMAP embeddings = global-z + per-trial baseline subtraction, backward windows, pacmap 0.8.0 (seed 0).\n');
fprintf(fid,'Input pulled from the shared windowed cache via getWindowed(''global'').\n');
fprintf(fid,'theta = angle (deg) between the Go and NoGo centroid-trajectory displacement vectors over t=0-0.8 s, per session.\n');
fprintf(fid,'n = %d sessions (Naive=%d, Trained=%d). All tests two-sided.\n\n', numel(th), sum(grp=="Naive"), sum(grp=="Trained"));
fprintf(fid,'Kruskal-Wallis:\n');
fprintf(fid,'  Group (Naive vs Trained), pooled over days : p = %.4g\n', p_kw_group);
fprintf(fid,'  Day/separation main effect                 : p = %.4g\n', p_kw_day);
fprintf(fid,'  Day effect within Naive                    : p = %.4g\n', p_kw_day_N);
fprintf(fid,'  Day effect within Trained                  : p = %.4g\n', p_kw_day_T);
fprintf(fid,'\nPer-day theta, median [IQR] deg (descriptive):\n');
for d=1:6, for g=["Naive","Trained"]
  v=th(day==d&grp==g); v=v(~isnan(v)); q=prctile(v,[25 50 75]);
  fprintf(fid,'  D%d %-7s %.1f [%.1f-%.1f] (n=%d)\n', d,g,q(2),q(1),q(3),numel(v));
end, end
fclose(fid); type(fullfile(statDir,'theta_stats.txt'));
fprintf('\nwrote theta deliverables: %s (CSVs), %s (stats)\n',dataDir,statDir);
