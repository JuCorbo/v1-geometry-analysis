%% sparsening_popmean_make.m  pop-mean panels + CSV + stats
% Reads ../cache/sparsening_popmean_pulled.mat (from sparsening_popmean_pull.m, shared cache).
% Timecourse tests: cluster-based permutation (groupSeparationTest, shared copy):
%   naive vs trained (Go, NoGo)  -> UNPAIRED (label shuffle)
%   Go vs NoGo (naive, trained)  -> PAIRED   (per-session sign-flip)
% all two-sided, |Hedges g|>=0.2, minSustained 2, 10000 perms.
codeDir=fileparts(mfilename('fullpath')); aDir=fileparts(codeDir);
figDir=fileparts(fileparts(aDir)); addpath(fullfile(figDir,'_sharedWindowedData','code'));
cacheDir=fullfile(aDir,'cache'); dataDir=fullfile(aDir,'data'); panDir=fullfile(aDir,'panels','main'); statDir=fullfile(aDir,'stats');
for d={dataDir,panDir,statDir}, if ~exist(d{1},'dir'), mkdir(d{1}); end, end
S=load(fullfile(cacheDir,'sparsening_popmean_pulled.mat'));
t=S.tvec(:)'; nT=numel(t);                     % PLOT/ANALYSE against window LAST frame (backward-window right edge)
cNaive=[0.85 0.65 0.13]; cTrain=[0.20 0.40 0.80];

%% PANEL C: rasters
clim=[-1.5 3]; m=256; zf=(0-clim(1))/(clim(2)-clim(1)); wi=max(2,min(m-1,round(zf*m)));
cmap=[[linspace(0,1,wi) linspace(1,1,m-wi)]',[linspace(0,1,wi) linspace(1,0,m-wi)]',[linspace(1,1,wi) linspace(1,0,m-wi)]'];
rasters={S.exampleNaive.matGo,S.exampleTrain.matGo;S.exampleNaive.matNo,S.exampleTrain.matNo};
rtitle={'Naive Go','Trained Go';'Naive NoGo','Trained NoGo'}; rfname={'naiveGo','trainedGo';'naiveNoGo','trainedNoGo'};
fR=figure('Color','w','Position',[80 80 700 640]); tlR=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
for i=1:2, for j=1:2
    ax=nexttile(tlR); M=rasters{i,j};
    imagesc(ax,t,1:size(M,1),M); set(ax,'YDir','reverse'); colormap(ax,cmap); caxis(ax,clim);
    xline(ax,0,'k-'); xline(ax,3,'k-'); xlim(ax,[-1 4]); title(ax,rtitle{i,j}); xlabel(ax,'Time (s)');
    if j==1, ylabel(ax,'Neuron'); end
    if j==2, cb=colorbar(ax); cb.Label.String='z-score'; end
    writetable(array2table(M,'VariableNames',compose('t_%.3f',t)), fullfile(dataDir,sprintf('fig1_e_raster_%s.csv',rfname{i,j})));
end, end
title(tlR,sprintf('Example sessions (D3, 30 deg): naive %s / trained %s',S.exampleNaive.name,S.exampleTrain.name),'Interpreter','none');
exportgraphics(fR,fullfile(panDir,'fig1_e_rasters.pdf'),'ContentType','vector');

%% PANEL D: example pop-mean timecourses
en=S.exampleNaive; et=S.exampleTrain;
band=@(ax,y,se,c,a) fill(ax,[t fliplr(t)],[y+se fliplr(y-se)],c,'FaceAlpha',a,'EdgeColor','none');
fD=figure('Color','w','Position',[80 80 1000 380]); tlD=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
axN=nexttile(tlD); hold(axN,'on'); box(axN,'off');
band(axN,en.tcGo,en.semGo,cNaive,0.18); band(axN,en.tcNo,en.semNo,cNaive,0.10);
plot(axN,t,en.tcGo,'-','Color',cNaive,'LineWidth',2.5); plot(axN,t,en.tcNo,'--','Color',cNaive,'LineWidth',2.5);
xline(axN,0,'k-'); xline(axN,3,'k-'); xlim(axN,[-1 4]); title(axN,'Naive (example)'); xlabel(axN,'Time (s)'); ylabel(axN,'Pop mean (z, win-int)');
axT=nexttile(tlD); hold(axT,'on'); box(axT,'off');
band(axT,et.tcGo,et.semGo,cTrain,0.18); band(axT,et.tcNo,et.semNo,cTrain,0.10);
plot(axT,t,et.tcGo,'-','Color',cTrain,'LineWidth',2.5); plot(axT,t,et.tcNo,'--','Color',cTrain,'LineWidth',2.5);
xline(axT,0,'k-'); xline(axT,3,'k-'); xlim(axT,[-1 4]); title(axT,'Trained (example)'); xlabel(axT,'Time (s)');
legend(axT,{'Go SEM','NoGo SEM','Go','NoGo'},'Box','off','Location','northeast'); linkaxes([axN axT],'y');
exportgraphics(fD,fullfile(panDir,'fig1_f_examplePopMeanTC.pdf'),'ContentType','vector');
writetable(table(t',en.tcGo',en.semGo',en.tcNo',en.semNo',et.tcGo',et.semGo',et.tcNo',et.semNo', ...
  'VariableNames',{'time_s','naiveGo_mean','naiveGo_sem','naiveNoGo_mean','naiveNoGo_sem','trainedGo_mean','trainedGo_sem','trainedNoGo_mean','trainedNoGo_sem'}), ...
  fullfile(dataDir,'fig1_f_examplePopMeanTC.csv'));

%% PANEL E: grand average + cluster-based permutation tests
mn=@(X) mean(X,1,'omitnan'); se=@(X) std(X,0,1,'omitnan')./sqrt(sum(~isnan(X),1));
mNaGo=mn(S.allNaiveGo); mNaNo=mn(S.allNaiveNo); mTrGo=mn(S.allTrainGo); mTrNo=mn(S.allTrainNo);
sNaGo=se(S.allNaiveGo); sNaNo=se(S.allNaiveNo); sTrGo=se(S.allTrainGo); sTrNo=se(S.allTrainNo);
nNa=size(S.allNaiveGo,1); nTr=size(S.allTrainGo,1);
dtCentre=mean(S.tvec(:)-S.tvec_centre(:));   % last frame - centre = 1 frame; report onset at first-sig-bin CENTRE
gst=@(A,B,pr) groupSeparationTest(A,B,t,'paired',pr,'correction','cluster','tail','both', ...
    'nPerms',10000,'minSustained',2,'minHedgesG',0.2,'onsetOffsetSec',dtCentre,'verbose',false);
rGo = gst(S.allNaiveGo, S.allTrainGo, false);   % naive vs trained, Go   (meanDiff = trained-naive)
rNo = gst(S.allNaiveNo, S.allTrainNo, false);   % naive vs trained, NoGo
rNa = gst(S.allNaiveGo, S.allNaiveNo, true);    % Go vs NoGo, naive   (paired; meanDiff = NoGo-Go)
rTr = gst(S.allTrainGo, S.allTrainNo, true);    % Go vs NoGo, trained (paired)
postIdx=t>=0;   % groupSeparationTest tests t>=0; map masks/g back onto the full 26-window axis
[sGo,gGo]=toFull(rGo,postIdx,nT); [sNo,gNo]=toFull(rNo,postIdx,nT);
[sNa,gNa]=toFull(rNa,postIdx,nT); [sTr,gTr]=toFull(rTr,postIdx,nT);

fE=figure('Color','w','Position',[80 80 900 460]); axE=axes(fE); hold(axE,'on'); box(axE,'off');
band(axE,mNaGo,sNaGo,cNaive,0.18); plot(axE,t,mNaGo,'-','Color',cNaive,'LineWidth',2.5);
band(axE,mNaNo,sNaNo,cNaive,0.10); plot(axE,t,mNaNo,'--','Color',cNaive,'LineWidth',2.5);
band(axE,mTrGo,sTrGo,cTrain,0.18); plot(axE,t,mTrGo,'-','Color',cTrain,'LineWidth',2.5);
band(axE,mTrNo,sTrNo,cTrain,0.10); plot(axE,t,mTrNo,'--','Color',cTrain,'LineWidth',2.5);
yl=ylim(axE); y0=yl(2);
drawSigRow(axE,t,sGo,y0*0.99,[0 0 0],'-',3);   % naive-vs-trained Go
drawSigRow(axE,t,sNo,y0*0.93,[0 0 0],'--',3);  % naive-vs-trained NoGo
xline(axE,0,'k-'); xline(axE,3,'k-'); xlim(axE,[-1 4]); xlabel(axE,'Time (s)'); ylabel(axE,'Pop mean (z, win-int)');
title(axE,sprintf('Grand average - Trained (n=%d), Naive (n=%d)',nTr,nNa));
exportgraphics(fE,fullfile(panDir,'fig1_g_grandAverage.pdf'),'ContentType','vector');
writetable(table(t',mNaGo',sNaGo',mNaNo',sNaNo',mTrGo',sTrGo',mTrNo',sTrNo', ...
  sGo(:),gGo(:), sNo(:),gNo(:), sNa(:),gNa(:), sTr(:),gTr(:), ...
  'VariableNames',{'time_s','naiveGo_mean','naiveGo_sem','naiveNoGo_mean','naiveNoGo_sem','trainedGo_mean','trainedGo_sem','trainedNoGo_mean','trainedNoGo_sem', ...
    'sig_naiveVsTrained_Go','g_naiveVsTrained_Go','sig_naiveVsTrained_NoGo','g_naiveVsTrained_NoGo', ...
    'sig_GoVsNoGo_naive','g_GoVsNoGo_naive','sig_GoVsNoGo_trained','g_GoVsNoGo_trained'}), ...
  fullfile(dataDir,'fig1_g_grandAverage_meanSEM.csv'));

%% stats section
stimIdx=t>=0&t<=3; susIdx=t>=1&t<=3;
fid=fopen(fullfile(statDir,'_popmean_stats.txt'),'w');
fprintf(fid,'--- Population mean activity (grand average across sessions) ---\n');
fprintf(fid,'Windowing: 3-frame backward, step 3, baseline z-score (sigma floor 0.12); plotted vs tvec (window last frame).\n');
fprintf(fid,'n sessions: naive Go=%d NoGo=%d | trained Go=%d NoGo=%d\n',size(S.allNaiveGo,1),size(S.allNaiveNo,1),size(S.allTrainGo,1),size(S.allTrainNo,1));
fprintf(fid,'Sustained mean +/- SEM over 1-3 s:\n');
fprintf(fid,'  naive Go   %.3f +/- %.3f | naive NoGo   %.3f +/- %.3f\n',mean(mNaGo(susIdx)),mean(sNaGo(susIdx)),mean(mNaNo(susIdx)),mean(sNaNo(susIdx)));
fprintf(fid,'  trained Go %.3f +/- %.3f | trained NoGo %.3f +/- %.3f\n',mean(mTrGo(susIdx)),mean(sTrGo(susIdx)),mean(mTrNo(susIdx)),mean(sTrNo(susIdx)));
gGo=hedgesg(mean(S.allNaiveGo(:,susIdx),2), mean(S.allTrainGo(:,susIdx),2));
gNo=hedgesg(mean(S.allNaiveNo(:,susIdx),2), mean(S.allTrainNo(:,susIdx),2));
fprintf(fid,'  Hedges g_av [BETWEEN-GROUP naive vs trained, UNPAIRED], sustained 1-3s: Go %+.2f | NoGo %+.2f  (g>0 => naive>trained)\n',gGo,gNo);
amps={mean(S.allNaiveGo(:,stimIdx),2),mean(S.allTrainGo(:,stimIdx),2),'Go'; mean(S.allNaiveNo(:,stimIdx),2),mean(S.allTrainNo(:,stimIdx),2),'NoGo'};
fprintf(fid,'\nStimulus-averaged population amplitude (per-session mean over 0-3 s):\n');
fprintf(fid,'  [BETWEEN-GROUP naive vs trained, UNPAIRED rank-sum; Hedges g_av; g>0 => naive>trained]\n');
for r=1:2
    aN=amps{r,1}; aN=aN(~isnan(aN)); aT=amps{r,2}; aT=aT(~isnan(aT));
    qN=prctile(aN,[25 50 75]); qT=prctile(aT,[25 50 75]);
    fprintf(fid,'  %-4s naive med %.3f[%.3f-%.3f] mean %.3f+/-%.3f SD (SEM %.3f) | trained med %.3f[%.3f-%.3f] mean %.3f+/-%.3f SD (SEM %.3f) | rank-sum p=%.3g Hedges g_av=%+.2f\n', ...
        amps{r,3}, qN(2),qN(1),qN(3), mean(aN),std(aN),std(aN)/sqrt(numel(aN)), ...
        qT(2),qT(1),qT(3), mean(aT),std(aT),std(aT)/sqrt(numel(aT)), ranksum(aN,aT), hedgesg(aN,aT));
end
fprintf(fid,'(n = %d naive, %d trained sessions.)\n', size(S.allNaiveGo,1), size(S.allTrainGo,1));
fprintf(fid,'\nWithin-group Go vs NoGo AMPLITUDE (per-session mean; PAIRED signed-rank; Hedges g_av; g>0 => NoGo>Go):\n');
wgrp={S.allNaiveGo,S.allNaiveNo,'naive'; S.allTrainGo,S.allTrainNo,'trained'};
wwin={susIdx,'1-3s(late)'; stimIdx,'0-3s(full)'};
for gi=1:2
    for wi=1:2
        mm=wwin{wi,1}; go=mean(wgrp{gi,1}(:,mm),2); no=mean(wgrp{gi,2}(:,mm),2);
        ok=~isnan(go)&~isnan(no); go=go(ok); no=no(ok); dd=no-go;
        fprintf(fid,'  %-7s %-11s Go med %.3f[%.3f-%.3f] mean %.3f+/-%.3f | NoGo med %.3f[%.3f-%.3f] mean %.3f+/-%.3f | diff(NoGo-Go) med %.3f mean %.3f+/-%.3f | signed-rank p=%.3g Hedges g_av=%+.2f\n', ...
            wgrp{gi,3}, wwin{wi,2}, median(go),prctile(go,25),prctile(go,75),mean(go),std(go), ...
            median(no),prctile(no,25),prctile(no,75),mean(no),std(no), median(dd),mean(dd),std(dd), signrank(go,no), hedgesg(no,go));
    end
end
fprintf(fid,'\nTimecourse: cluster-based permutation test (Maris & Oostenveld 2007), two-sided, 10000 perms,\n');
fprintf(fid,'|Hedges g|>=0.2, minSustained 2. Naive-vs-trained = UNPAIRED (label shuffle); Go-vs-NoGo = PAIRED (sign-flip).\n');
fprintf(fid,'Onset = CENTRE of the first sustained-significant window (s); timecourse plotted vs last-frame time. clusterP = fraction of perms with max cluster mass >= observed (floored 1e-4).\n');
cmp={'Naive vs Trained Go (unpaired)',rGo,sGo; 'Naive vs Trained NoGo (unpaired)',rNo,sNo; 'Go vs NoGo naive (paired)',rNa,sNa; 'Go vs NoGo trained (paired)',rTr,sTr};
for i=1:4
    r=cmp{i,2}; ss=cmp{i,3};
    fprintf(fid,'  %-34s onset t=%s s; clusterP=%.4g; sig %.0f%% of 0-3s bins\n', ...
        cmp{i,1}, num2str(r.onsetTime,'%.2f'), r.clusterPmin, 100*mean(ss(stimIdx)));
end
fprintf(fid,'\nExample sessions (D3, 30 deg): naive %s (%d neurons), trained %s (%d neurons)\n',S.exampleNaive.name,size(S.exampleNaive.matGo,1),S.exampleTrain.name,size(S.exampleTrain.matGo,1));
fclose(fid);
fprintf('sparsening_popmean_make done -> %s\n',aDir);

%% helpers
function g=hedgesg(x,y)
    x=x(~isnan(x)); y=y(~isnan(y)); n1=numel(x); n2=numel(y);
    sp=sqrt(((n1-1)*var(x)+(n2-1)*var(y))/(n1+n2-2));
    g=((mean(x)-mean(y))/sp)*(1-3/(4*(n1+n2)-9));
end
function [sig,g]=toFull(r,postIdx,nT)
    sig=false(1,nT); sig(postIdx)=r.sigMaskSust;
    g=nan(1,nT);     g(postIdx)=r.hedgesG;
end
function drawSigRow(ax,tv,sigMask,y,c,ls,lw)
    sigMask=sigMask(:)'; if ~any(sigMask), return; end
    tv=tv(:)'; d=diff([false sigMask false]); s=find(d==1); e=find(d==-1)-1;
    for i=1:numel(s), line(ax,[tv(s(i)) tv(e(i))],[y y],'Color',c,'LineStyle',ls,'LineWidth',lw); end
end
