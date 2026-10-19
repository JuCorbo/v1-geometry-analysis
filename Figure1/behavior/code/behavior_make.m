%% behavior_make.m — Figure 1 behavior: plot data (CSV) + stats from the cached intermediate
% Reads ../cache/behaviorSummary_data.mat (imaging D1-D6 hit/FA/d' for all trained
% animals) and ../data/fig1_c_training_exampleF07.csv (staged by makeF07TrainingCSV.m).
codeDir = fileparts(mfilename('fullpath'));
aDir    = fileparts(codeDir);
cacheDir= fullfile(aDir,'cache'); dataDir=fullfile(aDir,'data'); statDir=fullfile(aDir,'stats');

clip = @(p) max(min(p,.99),.01);

B = load(fullfile(cacheDir,'behaviorSummary_data.mat'));  % pHitAll,pFAAll,dpAll,animals
gngAng = [90 45 30 25 20 15];
[nDay,nA] = size(B.pHitAll);
animals = string(B.animals(:));

%% CSV 1: imaging D1-D6 per animal
rows = {};
for d=1:nDay, for a=1:nA
    rows(end+1,:)={d,gngAng(d),char(animals(a)),B.pHitAll(d,a),B.pFAAll(d,a),B.dpAll(d,a)}; %#ok<SAGROW>
end, end
writetable(cell2table(rows,'VariableNames',{'day','gngAngle_deg','animal','hitRate','faRate','dprime'}), ...
    fullfile(dataDir,'fig1_d_ROC_perAnimal_D1D6.csv'));

%% CSV 2: group mean +/- SEM per day
mHit=mean(B.pHitAll,2,'omitnan'); mFA=mean(B.pFAAll,2,'omitnan'); mDp=mean(B.dpAll,2,'omitnan');
nDp=sum(~isnan(B.dpAll),2);
sHit=std(B.pHitAll,0,2,'omitnan')./sqrt(sum(~isnan(B.pHitAll),2));
sFA =std(B.pFAAll ,0,2,'omitnan')./sqrt(sum(~isnan(B.pFAAll),2));
sDp =std(B.dpAll  ,0,2,'omitnan')./sqrt(nDp);
writetable(table((1:nDay)',gngAng',mHit,sHit,mFA,sFA,mDp,sDp,nDp, ...
  'VariableNames',{'day','gngAngle_deg','meanHit','semHit','meanFA','semFA','meanDprime','semDprime','nAnimals'}), ...
  fullfile(dataDir,'fig1_d_ROC_groupMeanSEM_D1D6.csv'));

%% STATS
nCeil=128; pGoCeil=mean(B.pHitAll(:),'omitnan');
CEIL=norminv(clip(binoinv(0.775,nCeil,pGoCeil)/nCeil))-norminv(clip(binoinv(0.225,nCeil,0.5)/nCeil));
sep=repmat(gngAng(:),1,nA); okv=~isnan(B.dpAll);
[rhoSep,pSep]=corr(sep(okv),B.dpAll(okv),'type','Spearman');
nAbove=sum(mDp>CEIL);
fid=fopen(fullfile(statDir,'behavior_stats.txt'),'w');
fprintf(fid,'FIGURE 1 - BEHAVIOR: statistics\nGenerated %s\n',datestr(now));
fprintf(fid,'Source (cache): behaviorSummary_data.mat; fig1_c_training_exampleF07.csv\n\n');
fprintf(fid,'Trained cohort during imaging: n = %d animals (%s)\n\n',nA,strjoin(cellstr(animals),', '));
fprintf(fid,'Per-day performance (group mean +/- SEM across animals):\n');
fprintf(fid,'%-4s %-8s %-16s %-16s %-16s %-4s\n','day','sep(deg)','hitRate','faRate','dprime','n');
for d=1:nDay
  fprintf(fid,'D%-3d %-8d %.3f +/- %.3f   %.3f +/- %.3f   %.3f +/- %.3f   %d\n', ...
    d,gngAng(d),mHit(d),sHit(d),mFA(d),sFA(d),mDp(d),sDp(d),nDp(d));
end
fprintf(fid,'\nChance ceiling (best random d''): %.3f\n',CEIL);
fprintf(fid,'  method: favourable 77.5%%/22.5%% percentile corner, n=%d, pGo=%.3f (Corbo 2025).\n',nCeil,pGoCeil);
fprintf(fid,'  day-means exceeding ceiling: %d / %d\n',nAbove,nDay);
fprintf(fid,'\nDifficulty effect (Spearman d'' vs separation): rho=%.3f, p=%.3g, n=%d\n',rhoSep,pSep,sum(okv(:)));
E=readtable(fullfile(dataDir,'fig1_c_training_exampleF07.csv'));
fprintf(fid,'\nExample training animal F07: %d sessions, session d'' %.2f -> %.2f\n', ...
  height(E),E.dprime(1),E.dprime(end));
fclose(fid);
fprintf('behavior_make done -> %s\n',aDir);
