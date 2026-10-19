%% fractionActive_make.m  active-neuron-fraction panel (Go+NoGo combined)
% Fig 1h: fraction of activated neurons for ALL task trials (Go+NoGo pooled),
% naive vs trained, in three 1 s windows. Reads ../cache/fractionActive_windowed.mat
% (from fractionActive_pull.m, shared cache). Also merges pop-mean + fraction stats into
% sparsening_stats.txt.
codeDir=fileparts(mfilename('fullpath')); aDir=fileparts(codeDir);
cacheDir=fullfile(aDir,'cache'); dataDir=fullfile(aDir,'data'); panDir=fullfile(aDir,'panels','main'); statDir=fullfile(aDir,'stats');
L=load(fullfile(cacheDir,'fractionActive_windowed.mat'));   % fracAct_binom,sessName_binom,gngAng,colNames,condLab,NR
frac=L.fracAct_binom; nm=L.sessName_binom; gngAng=L.gngAng; colNames=L.colNames;
nDay=size(frac,2); grpName={'naive','trained'};
cNaive=[0.85 0.65 0.13]; cTrain=[0.20 0.40 0.80];

%% tidy per-session CSV (all 9 columns: All/Go/NoGo x 3 windows)
rows={};
for gi=1:2, for di=1:nDay
    M=frac{gi,di}; s=nm{gi,di};
    for r=1:size(M,1), rows(end+1,:)=[{grpName{gi},di,gngAng(di),s{r}},num2cell(M(r,:))]; end %#ok<SAGROW>
end, end
writetable(cell2table(rows,'VariableNames',[{'group','day','gngAngle_deg','session'},colNames]), ...
    fullfile(dataDir,'fig1_h_fractionActive_perSession.csv'));

%% pool across days
pool={[],[]}; for gi=1:2, for di=1:nDay, pool{gi}=[pool{gi};frac{gi,di}]; end, end
allCols=1:3; winLab={'0-1s','1-2s','2-3s'};
pMain=nan(1,3); for w=1:3, a=pool{1}(:,allCols(w)); a=a(~isnan(a)); b=pool{2}(:,allCols(w)); b=b(~isnan(b)); pMain(w)=ranksum(a,b); end
qMain=bh(pMain);

%% MAIN panel: combined fraction, naive vs trained, 3 windows
fF=figure('Color','w','Position',[80 80 380 380]); ax=axes(fF); hold(ax,'on'); box(ax,'off');
ymax=0;
for w=1:3
    a=pool{1}(:,allCols(w)); a=a(~isnan(a)); b=pool{2}(:,allCols(w)); b=b(~isnan(b));
    xa=w-0.18; xb=w+0.18; ymax=max([ymax;a;b]);
    scatter(ax,xa+0.05*randn(numel(a),1),a,9,cNaive,'filled','MarkerFaceAlpha',0.35);
    scatter(ax,xb+0.05*randn(numel(b),1),b,9,cTrain,'filled','MarkerFaceAlpha',0.35);
    errorbar(ax,xa,mean(a),std(a)/sqrt(numel(a)),'o','Color',cNaive,'MarkerFaceColor',cNaive,'LineWidth',1.6,'CapSize',9);
    errorbar(ax,xb,mean(b),std(b)/sqrt(numel(b)),'o','Color',cTrain,'MarkerFaceColor',cTrain,'LineWidth',1.6,'CapSize',9);
    s=stars(qMain(w)); if ~strcmp(s,'n.s.'), text(ax,w,max([a;b])*1.03,s,'HorizontalAlignment','center','FontSize',13); end
end
set(ax,'XTick',1:3,'XTickLabel',winLab); xlim(ax,[0.5 3.5]); ylim(ax,[0 ymax*1.12]);
xlabel(ax,'Stimulus window'); ylabel(ax,'Fraction of neurons activated');
title(ax,'Go+NoGo combined (naive gold / trained blue)');
exportgraphics(fF,fullfile(panDir,'fig1_h_fractionActive.pdf'),'ContentType','vector');

%% stats section (report All main + Go + NoGo)
fid=fopen(fullfile(statDir,'_fraction_active_stats.txt'),'w');
fprintf(fid,'--- Active-neuron fraction: MAIN = All (Go+NoGo pooled trials) ---\n');
fprintf(fid,'Comparison: BETWEEN-GROUP naive vs trained (UNPAIRED); Wilcoxon rank-sum per window; BH-FDR across 3 windows.\n');
fprintf(fid,'Effect size: Hedges g_av (pooled within-condition SD); g>0 => naive>trained.\n');
fprintf(fid,'Neuron activated if window P(APrE>=1) > baseline Clopper-Pearson 95%% upper CI (alpha=0.05).\n\n');
fprintf(fid,'%-11s  %-20s  %-20s  %-11s %-9s %-9s\n','window','naive med[IQR]','trained med[IQR]','ranksum p','BH q','Hedges g_av');
for c=1:3
    a=pool{1}(:,c); a=a(~isnan(a)); b=pool{2}(:,c); b=b(~isnan(b));
    qa=prctile(a,[25 50 75]); qb=prctile(b,[25 50 75]);
    fprintf(fid,'%-11s  %.3f[%.3f-%.3f]  %.3f[%.3f-%.3f]  %-11.3g %-9.3g %+.2f\n', ...
        colNames{c},qa(2),qa(1),qa(3),qb(2),qb(1),qb(3),pMain(c),qMain(c),hedgesg(a,b));
end
fprintf(fid,'(Hedges g_av > 0 => naive > trained; bands: <0.2 negligible, 0.2-0.4 small, 0.4-0.6 medium, >0.6 large.)\n');
fprintf(fid,'\n(mean +/- SD):\n');
for c=1:3
    a=pool{1}(:,c); a=a(~isnan(a)); b=pool{2}(:,c); b=b(~isnan(b));
    fprintf(fid,'  %-11s naive %.3f +/- %.3f | trained %.3f +/- %.3f\n',colNames{c},mean(a),std(a),mean(b),std(b));
end
fprintf(fid,'\nStimulus-averaged responsive fraction (per-session mean over 0-1/1-2/2-3 s):\n');
fprintf(fid,'  [BETWEEN-GROUP naive vs trained, UNPAIRED rank-sum; Hedges g_av; g>0 => naive>trained]\n');
segNames={'All','Go','NoGo'}; segCols={1:3,4:6,7:9};
for s=1:3
    aN=mean(pool{1}(:,segCols{s}),2); aN=aN(~isnan(aN));
    aT=mean(pool{2}(:,segCols{s}),2); aT=aT(~isnan(aT));
    qN=prctile(aN,[25 50 75]); qT=prctile(aT,[25 50 75]);
    fprintf(fid,'  %-5s naive med %.3f[%.3f-%.3f] mean %.3f+/-%.3f SD (SEM %.3f) | trained med %.3f[%.3f-%.3f] mean %.3f+/-%.3f SD (SEM %.3f) | rank-sum p=%.3g Hedges g_av=%+.2f\n', ...
        segNames{s}, qN(2),qN(1),qN(3), mean(aN),std(aN),std(aN)/sqrt(numel(aN)), ...
        qT(2),qT(1),qT(3), mean(aT),std(aT),std(aT)/sqrt(numel(aT)), ranksum(aN,aT), hedgesg(aN,aT));
end
fprintf(fid,'(n = %d naive, %d trained sessions.)\n', size(pool{1},1), size(pool{2},1));
fclose(fid);

%% merge into sparsening_stats.txt
outF=fullfile(statDir,'sparsening_stats.txt'); fo=fopen(outF,'w');
fprintf(fo,'FIGURE 1 - SPARSENING: statistics\nGenerated %s\n\n',datestr(now));
fprintf(fo,'CONVENTIONS\n');
fprintf(fo,'  Signal = baseline z-scored (sigma floor 0.12), window-integrated APrEs; 3-frame backward windows, plotted vs tvec (window last frame).\n');
fprintf(fo,'  Timecourse (pop-mean): cluster-based permutation (Maris & Oostenveld). Between-group naive-vs-trained = UNPAIRED; within-group Go-vs-NoGo = PAIRED (sign-flip).\n');
fprintf(fo,'  Windowed (fraction active, amplitudes): rank-sum (between-group), signed-rank (within-group), BH-FDR.\n');
fprintf(fo,'  Effect size = Hedges g_av; bands: <0.2 negligible, 0.2-0.4 small, 0.4-0.6 medium, >0.6 large. Median [IQR]; means +/- SD.\n\n');
fprintf(fo,'========== (A) POPULATION-MEAN AMPLITUDE ==========\n');
fp=fullfile(statDir,'_popmean_stats.txt'); if exist(fp,'file'), fprintf(fo,'%s\n\n',fileread(fp)); end
fprintf(fo,'========== (B) FRACTION OF RESPONSIVE NEURONS ==========\n');
fp=fullfile(statDir,'_fraction_active_stats.txt'); if exist(fp,'file'), fprintf(fo,'%s\n\n',fileread(fp)); end
fclose(fo);
fprintf('fractionActive_make (MAIN=All) done; merged -> %s\n',outF);

%% helpers
function g=hedgesg(x,y)
    x=x(~isnan(x)); y=y(~isnan(y)); n1=numel(x); n2=numel(y);
    sp=sqrt(((n1-1)*var(x)+(n2-1)*var(y))/(n1+n2-2));
    g=((mean(x)-mean(y))/sp)*(1-3/(4*(n1+n2)-9));
end
function s=stars(q), if q<0.001,s='***';elseif q<0.01,s='**';elseif q<0.05,s='*';else,s='n.s.';end, end
function qv=bh(pv)
    pv=pv(:); m=numel(pv); [ps,ix]=sort(pv); adj=ps.*m./(1:m)';
    adj=min(1,flipud(cummin(flipud(adj)))); qv=zeros(m,1); qv(ix)=adj; qv=reshape(qv,1,[]);
end
