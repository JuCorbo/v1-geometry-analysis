%% Whole-database sigma diagnostics: pooled sigma / z-noise, trained vs naive.
codeDir    = fileparts(mfilename('fullpath'));
figRoot    = fileparts(fileparts(fileparts(codeDir)));
sharedCode = fullfile(figRoot,'_sharedWindowedData','code');
sharedData = fullfile(figRoot,'_sharedWindowedData','data');
addpath(sharedCode);
files = dir(fullfile(sharedData,'*','D*','*_win.mat'));
Sig={}; Zsd={}; Maz={}; Grp={};
sFrac=[]; sGrp=strings(0,1); sN=[];
for k = 1:numel(files)
    S = load(fullfile(files(k).folder,files(k).name));
    [Xnf,info] = getWindowed(S,'baseline_nofloor');
    sig = double(S.sd_base_unfloored(:)); tv = info.tvec(:); post = tv>0;
    Zp = Xnf(:,:,post); nTr=size(Zp,1); nC=size(Zp,2); nP=sum(post);
    M = reshape(permute(Zp,[1 3 2]), nTr*nP, nC);
    Sig{end+1}=sig; Zsd{end+1}=std(M,0,1,'omitnan')'; Maz{end+1}=mean(abs(M),1,'omitnan')'; %#ok<*SAGROW>
    [p1,~]=fileparts(files(k).folder); [~,g]=fileparts(p1);
    Grp{end+1}=repmat(string(g),nC,1);
    sFrac(end+1)=mean(sig<0.12); sGrp(end+1)=string(g); sN(end+1)=numel(info.baselineBins)*nTr;
end
Sig=cat(1,Sig{:}); Zsd=cat(1,Zsd{:}); Maz=cat(1,Maz{:}); Grp=cat(1,Grp{:});
Zsd(~isfinite(Zsd))=NaN; Maz(~isfinite(Maz))=NaN;
isT=Grp=="trained"; isN=Grp=="naive";
cT=[0.20 0.40 0.80]; cN=[0.85 0.65 0.13]; floors=0.12;   % adopted floor

% --- explicit log-binned medians of a quantity vs sigma ---
Sv = Sig(isfinite(Sig)&Sig>0);   % prctile returns NaN if any NaN present -> filter first
edges = logspace(log10(prctile(Sv,0.3)), log10(prctile(Sv,99.7)), 22); nb=numel(edges)-1;
ctr = sqrt(edges(1:end-1).*edges(2:end));
mTz=nan(1,nb); mNz=nan(1,nb); q1Tz=nan(1,nb); q3Tz=nan(1,nb); q1Nz=nan(1,nb); q3Nz=nan(1,nb); mTm=nan(1,nb); mNm=nan(1,nb);
for i=1:nb
    inb = Sig>=edges(i) & Sig<edges(i+1);
    zt=Zsd(inb&isT); zn=Zsd(inb&isN); mt=Maz(inb&isT); mn=Maz(inb&isN);
    if any(~isnan(zt)), mTz(i)=median(zt,'omitnan'); q1Tz(i)=prctile(zt,25); q3Tz(i)=prctile(zt,75); end
    if any(~isnan(zn)), mNz(i)=median(zn,'omitnan'); q1Nz(i)=prctile(zn,25); q3Nz(i)=prctile(zn,75); end
    if any(~isnan(mt)), mTm(i)=median(mt,'omitnan'); end
    if any(~isnan(mn)), mNm(i)=median(mn,'omitnan'); end
end

fig=figure('Color','w','Position',[40 40 1650 1150]);
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
title(tl,sprintf('Whole-database \\sigma diagnostics  |  %d neurons, %d sessions  |  Trained (blue) vs Naive (gold)',numel(Sig),numel(files)),'FontWeight','bold');

% (a) pooled sigma distribution
ax=nexttile; hold(ax,'on'); se=logspace(log10(min(Sig(Sig>0))),log10(max(Sig)),50);
histogram(ax,Sig(isT),se,'Normalization','probability','FaceColor',cT,'FaceAlpha',0.45,'EdgeColor','none');
histogram(ax,Sig(isN),se,'Normalization','probability','FaceColor',cN,'FaceAlpha',0.45,'EdgeColor','none');
set(ax,'XScale','log'); xlim(ax,[se(1) se(end)]); for f=floors, xline(ax,f,'-','Color',[.3 .3 .3],'LineWidth',1.3,'Label',sprintf('%.2f',f)); end
xlabel(ax,'\sigma (unfloored baseline SD)'); ylabel(ax,'prob'); title(ax,'(a) pooled \sigma distribution');
text(ax,0.02,0.97,sprintf('%%<0.12: Tr %.0f%% Na %.0f%%',100*mean(Sig(isT)<0.12),100*mean(Sig(isN)<0.12)),'Units','normalized','VerticalAlignment','top','FontSize',9);

% (b) sigma vs single-trial z SD (binned median + IQR)
ax=nexttile; hold(ax,'on');
gd=~isnan(mTz); fill(ax,[ctr(gd) fliplr(ctr(gd))],[q1Tz(gd) fliplr(q3Tz(gd))],cT,'FaceAlpha',0.15,'EdgeColor','none');
gd=~isnan(mNz); fill(ax,[ctr(gd) fliplr(ctr(gd))],[q1Nz(gd) fliplr(q3Nz(gd))],cN,'FaceAlpha',0.15,'EdgeColor','none');
plot(ax,ctr,mTz,'-','Color',cT,'LineWidth',2.4); plot(ax,ctr,mNz,'-','Color',cN,'LineWidth',2.4);
ref=0.9*(0.2./ctr); plot(ax,ctr,ref,'k--','LineWidth',1); text(ax,ctr(4),ref(4),' \propto 1/\sigma','FontSize',9);
set(ax,'XScale','log','YScale','log'); xlim(ax,[edges(1) edges(end)]);
for f=floors, xline(ax,f,'-','Color',[.3 .3 .3],'LineWidth',1.3); end
xlabel(ax,'\sigma (log)'); ylabel(ax,'single-trial z SD (median, IQR band)'); title(ax,'(b) \sigma vs single-trial z SD  [the DR noise]');

% (c) sigma vs mean |z| (binned median)
ax=nexttile; hold(ax,'on');
plot(ax,ctr,mTm,'-','Color',cT,'LineWidth',2.4); plot(ax,ctr,mNm,'-','Color',cN,'LineWidth',2.4);
set(ax,'XScale','log','YScale','log'); xlim(ax,[edges(1) edges(end)]);
for f=floors, xline(ax,f,'-','Color',[.3 .3 .3],'LineWidth',1.3); end
xlabel(ax,'\sigma (log)'); ylabel(ax,'mean post-stim |z| (median)'); title(ax,'(c) \sigma vs mean |z|');

% (d) per-session spread of the <0.10 fraction
ax=nexttile; hold(ax,'on');
histogram(ax,100*sFrac(sGrp=="trained"),0:5:80,'Normalization','probability','FaceColor',cT,'FaceAlpha',0.45,'EdgeColor','none');
histogram(ax,100*sFrac(sGrp=="naive"),0:5:80,'Normalization','probability','FaceColor',cN,'FaceAlpha',0.45,'EdgeColor','none');
xlabel(ax,'per-session %% neurons with \sigma<0.12'); ylabel(ax,'prob'); title(ax,'(d) per-session spread of <0.12 fraction');
text(ax,0.98,0.97,sprintf('n\\_base %d-%d\nTr med %.0f%%, Na med %.0f%%',min(sN),max(sN),median(100*sFrac(sGrp=="trained")),median(100*sFrac(sGrp=="naive"))),'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','FontSize',9);

outDir=fullfile(fileparts(codeDir),'panels','supp');
if ~exist(outDir,'dir'), mkdir(outDir); end
outPng=fullfile(outDir,'sigma_floor_diagnostic.png');
exportgraphics(fig,outPng,'Resolution',150,'BackgroundColor','white');
exportgraphics(fig,strrep(outPng,'.png','.pdf'),'ContentType','vector','BackgroundColor','white');
fprintf('saved %s\n',outPng);
fprintf('binned zSD trained: '); fprintf('%.2f ',mTz); fprintf('\nctr: '); fprintf('%.3f ',ctr); fprintf('\n');
