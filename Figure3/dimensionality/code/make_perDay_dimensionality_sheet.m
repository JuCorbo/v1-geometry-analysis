%% make_perDay_dimensionality_sheet.m
%  Per-day dimensionality supplement: rows = measure (participation ratio D_PR,
%  k80, tail slope), columns = day D1..D6. Each cell = naive/trained per-session
%  traces (thin) + group means (thick), x-axis from baseline (~ -1 s), one shared
%  y-scale per row. Reads the Figure 3 dimensionality cache (dimMetrics), which
%  holds per-session/day/window D_PR/k80/tailSlope for Go and NoGo across 5
%  subsample rounds; each session timecourse = mean over Go/NoGo and rounds.
clear; clc;
codeDir=fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
anaDir=fileparts(codeDir);
D=load(fullfile(anaDir,'cache','dimMetrics_baseline012.mat'));
dm=D.dimMetrics; tvec=D.tvec(:)'; nW=numel(tvec);
suppDir=fullfile(anaDir,'panels','supp'); if ~exist(suppDir,'dir'), mkdir(suppDir); end
% dbi: 1 = naive (BatuhanADNaive), 2 = trained (BatuhanAngDist)
colNaive=[0.85 0.65 0.13]; colTrain=[0.20 0.40 0.80]; thinAlpha=0.22; thinLW=0.8; thickLW=2.4;
dayLab=arrayfun(@(x)sprintf('D%d',x),1:6,'uni',0);

measures={'D_PR','Participation ratio (D_{PR})'; 'k80','k80'; 'tailSlope','Tail slope'};
nM=size(measures,1); nD=6;

fig=figure('Color','w','Position',[40 40 1700 950]);
tl=tiledlayout(nM,nD,'TileSpacing','compact','Padding','compact');
sgtitle(tl,'Neural dimensionality, per day (rows = measure, columns = day) | naive vs trained','FontWeight','bold');

for r=1:nM
  mf=measures{r,1};
  TC=cell(2,nD);
  for dbi=1:2, for d=1:nD, TC{dbi,d}=sessTC(dm,dbi,d,mf,nW); end, end
  % per-row shared y-limits over all days/both groups
  allY=[]; for dbi=1:2, for d=1:nD, allY=[allY; TC{dbi,d}(:)]; end, end %#ok<AGROW>
  allY=allY(isfinite(allY)); yl=[min(allY) max(allY)];
  if isempty(allY)||yl(1)==yl(2), yl=[0 1]; end
  pad=0.05*(yl(2)-yl(1)); yl=[yl(1)-pad yl(2)+pad];

  for d=1:nD
    ax=nexttile(tl,(r-1)*nD+d); hold(ax,'on'); box(ax,'off');
    plotSessions(ax,tvec,TC{1,d},colNaive,thinAlpha,thinLW);
    plotSessions(ax,tvec,TC{2,d},colTrain,thinAlpha,thinLW);
    if ~isempty(TC{1,d}), plot(ax,tvec,mean(TC{1,d},1,'omitnan'),'-','Color',colNaive,'LineWidth',thickLW); end
    if ~isempty(TC{2,d}), plot(ax,tvec,mean(TC{2,d},1,'omitnan'),'-','Color',colTrain,'LineWidth',thickLW); end
    xline(ax,0,':','Color',[0.6 0.6 0.6]); xlim(ax,[tvec(1) tvec(end)]); ylim(ax,yl); set(ax,'FontSize',7);
    if r==1, title(ax,dayLab{d},'Interpreter','none'); end
    if d==1, ylabel(ax,measures{r,2}); else, set(ax,'YTickLabel',{}); end
    if r<nM, set(ax,'XTickLabel',{}); else, xlabel(ax,'Time (s)'); end
    if r==1 && d==nD
      h1=plot(ax,nan,nan,'-','Color',blendWhite(colNaive,1-thinAlpha),'LineWidth',thinLW);
      h2=plot(ax,nan,nan,'-','Color',blendWhite(colTrain,1-thinAlpha),'LineWidth',thinLW);
      h3=plot(ax,nan,nan,'-','Color',colNaive,'LineWidth',thickLW);
      h4=plot(ax,nan,nan,'-','Color',colTrain,'LineWidth',thickLW);
      legend(ax,[h1 h2 h3 h4],{'naive sessions','trained sessions','naive mean','trained mean'},'FontSize',6,'Box','off','Location','northeast');
    end
  end
end

base=fullfile(suppDir,'supp_dimensionality_perDay');
exportgraphics(fig,[base '.pdf'],'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,[base '.png'],'Resolution',160,'BackgroundColor','white');
savefig(fig,[base '.fig']);
fprintf('Saved per-day dimensionality sheet to %s.(pdf/png/fig)\n',base);

%% helpers
function TC=sessTC(dm,dbi,day,mf,nW)
  TC=[]; nSlots=size(dm,3); nRep=size(dm,4);
  for s=1:nSlots
    reps=[];
    for r=1:nRep
      e=dm(dbi,day,s,r); if isempty(e.go), continue; end
      g=arrayfun(@(x)x.(mf),e.go(:))'; n=arrayfun(@(x)x.(mf),e.nogo(:))';
      if numel(g)~=nW||numel(n)~=nW, continue; end
      reps(end+1,:)=0.5*(g+n); %#ok<AGROW>
    end
    if ~isempty(reps), TC(end+1,:)=mean(reps,1,'omitnan'); end %#ok<AGROW>
  end
end
function plotSessions(ax,tvec,M,baseCol,alphaLike,lw)
  if isempty(M), return; end
  c=blendWhite(baseCol,1-alphaLike);
  for i=1:size(M,1), plot(ax,tvec,M(i,:),'-','Color',c,'LineWidth',lw); end
end
function c=blendWhite(baseCol,fracWhite)
  fracWhite=min(max(fracWhite,0),1); c=(1-fracWhite)*baseCol+fracWhite*[1 1 1];
end
