%% make_windowed_cache.m -- build the shared windowed intermediate from a Zenodo-format CSV export
%
% CONVENTIONS:
%   backward windows, every bin a full 3 frames : bins = (c-2):c for c = 9:3:84
%                                                 -> {7,8,9} {10,11,12} ... {82,83,84}
%   stored sd_base is floored at 0.05; analyses use getWindowed, which applies
%     the paper-wide 0.12 floor to the stored sd_base_unfloored instead
%   raw (un-normalised) sums stored; normalisation constants stored alongside
%
% Usage:  GROUP='trained'; IN='ds1_filtered.csv'; make_windowed_cache
%
% Input CSV columns: Experiment,Trial,Cell,Behav_Cond,Visual_Stim,TraceSpikes(hex)
% (already QC- and task-block-filtered by the stage-1 awk pass)

if ~exist('GROUP','var'), GROUP = 'trained'; end
% SCR = staging folder holding the stage-1 awk extraction of the Zenodo deposit
% (ds1_filtered.csv / nv_ordered.csv etc. -- see header). Not distributed: the
% windowed cache in ../data is the product of this script. Set SCR to rebuild.
if ~exist('SCR','var'), error('Set SCR to the stage-1 CSV extraction folder (see header comments).'); end
if ~exist('OUTROOT','var'), OUTROOT = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data'); end
if ~exist(OUTROOT,'dir'), mkdir(OUTROOT); end

% group-specific column layout -- the naive export has an extra `layer` column and
% only the *_spikes_1 tuning set, so indices differ from trained.
switch GROUP
    case 'trained'
        DBNAME='BatuhanAngDist'; PREFCOL=27; OSICOL=30;              % *_spikes_2
        IN_CSV  = fullfile(SCR,'ds1_filtered.csv');
        CELLTBL = fullfile(SCR,'ds1','DATASET1_CellTable.csv');
        TRIALTBL= fullfile(SCR,'ds1_trials_small.csv');
    case 'naive'
        DBNAME='BatuhanADNaive'; PREFCOL=10; OSICOL=13;              % *_spikes_1
        IN_CSV  = fullfile(SCR,'nv_ordered.csv');
        CELLTBL = fullfile(SCR,'nv','Naive_CellTable.csv');
        TRIALTBL= fullfile(SCR,'nv_trials_small.csv');
    otherwise, error('GROUP must be trained or naive');
end
fprintf('GROUP=%s  DB=%s  prefCol=%d osiCol=%d\n  in: %s\n', GROUP, DBNAME, PREFCOL, OSICOL, IN_CSV);

%% ---------------- window definition ----------------
dt=0.0646; rasterLen=140; winFrames=3; stepFrames=3;
xval=linspace(0,(rasterLen-1)*dt,rasterLen)-1.5;
tcWindowSec=[-1 4];
tcIdx=find(xval>=tcWindowSec(1) & xval<=tcWindowSec(2));
startIdx=tcIdx(1); endIdx=tcIdx(end);
centers=startIdx:stepFrames:endIdx;
bins=arrayfun(@(c) (c-winFrames+1):c, centers, 'uni', 0);   % ALL full width
nWin=numel(bins);
tvec=xval(centers);                          % last-frame anchor (matches the scripts)
tvec_centre=xval(centers-(winFrames-1)/2);   % TRUE centre of a backward window
baselineBins=find(tvec<0);
baselineWindowSec=[-1.5 0];
baselineWindow=find(xval>=baselineWindowSec(1) & xval<=baselineWindowSec(2));
SD_FLOOR=0.05; GLOB_FLOOR=1e-8;
assert(all(cellfun(@numel,bins)==winFrames), 'all bins must be full width');
fprintf('bins: %d x %d frames, %s .. %s ; baseline bins %s (frames %d..%d)\n', ...
    nWin, winFrames, mat2str(bins{1}), mat2str(bins{end}), mat2str(baselineBins), ...
    bins{baselineBins(1)}(1), bins{baselineBins(end)}(end));
fprintf('baselineWindow for subtraction: frames %d..%d\n', baselineWindow(1), baselineWindow(end));

%% ---------------- per-cell OSI / pref-orientation ----------------
fprintf('reading CellTable ...\n');
fid=fopen(CELLTBL,'r'); fgetl(fid);
cellKey={}; cellOSI=[]; cellPref=[];
while true
    L=fgetl(fid); if ~ischar(L), break; end
    p=strsplit(L,',');
    if numel(p)<max(PREFCOL,OSICOL), continue; end
    cellKey{end+1,1}=[p{2} '|' p{1}]; %#ok<AGROW>
    cellPref(end+1,1)=str2double(p{PREFCOL}); %#ok<AGROW>
    cellOSI(end+1,1)=str2double(p{OSICOL});   %#ok<AGROW>
end
fclose(fid);
cellMap=containers.Map(cellKey,(1:numel(cellKey))');
fprintf('  %d cells indexed\n', numel(cellKey));

%% ---------------- per-trial behaviour ----------------
fprintf('reading TrialTable ...\n');
T2=readtable(TRIALTBL,'ReadVariableNames',false,'Format','%f%s%s%s%f%f%s','Delimiter',',');
T2.Properties.VariableNames={'Trial','Experiment','Behav_Cond','Block','Visual_Stim','Response','Outcome'};
trialMap=containers.Map(strcat(T2.Experiment,'|',string(T2.Trial)), (1:height(T2))');

%% ---------------- stream ----------------
man=cell(0,8);
fid=fopen(IN_CSV,'r'); assert(fid>0,'cannot open %s',IN_CSV);
curExp=''; curDay=''; bT=[];bC=[];bV=[];bH={}; nS=0; t0=tic;
seenExp=containers.Map();   % guard: the reader assumes each session is one contiguous block.
                            % If the input is not sorted by Experiment a session would be written
                            % twice and silently truncated -- fail loudly instead.
while true
    L=fgetl(fid); atEOF=~ischar(L);
    if ~atEOF
        k=find(L==',',5);
        ex=L(1:k(1)-1); tri=str2double(L(k(1)+1:k(2)-1)); cel=str2double(L(k(2)+1:k(3)-1));
        day=L(k(3)+1:k(4)-1); vst=str2double(L(k(4)+1:k(5)-1)); hx=L(k(5)+1:end);
    end
    if atEOF || (~isempty(curExp) && ~strcmp(ex,curExp))
        r=writeSession(curExp,curDay,bT,bC,bV,bH,DBNAME,GROUP,OUTROOT, ...
            rasterLen,endIdx,bins,nWin,baselineBins,baselineWindow,SD_FLOOR,GLOB_FLOOR, ...
            centers,tvec,tvec_centre,dt,winFrames,stepFrames,tcWindowSec, ...
            cellMap,cellOSI,cellPref,trialMap,T2);
        if ~isempty(r)
            nS=nS+1; man(end+1,:)=r; %#ok<AGROW>
            if mod(nS,15)==0||nS==1, fprintf('[%3d] %-14s %s  %dx%d  (%.0f s)\n',nS,r{2},r{3},r{4},r{5},toc(t0)); end
        end
        bT=[];bC=[];bV=[];bH={};
    end
    if atEOF, break; end
    if ~strcmp(ex,curExp)
        if isKey(seenExp,ex)
            error('make_windowed_cache:notContiguous', ...
                ['session %s appears in more than one block of %s. The reader requires each ' ...
                 'session to be contiguous; sort the input by Experiment first.'], ex, IN_CSV);
        end
        seenExp(ex)=true;
    end
    curExp=ex; curDay=day;
    bT(end+1,1)=tri; bC(end+1,1)=cel; bV(end+1,1)=vst; bH{end+1,1}=hx; %#ok<AGROW>
end
fclose(fid);

M=cell2table(man,'VariableNames',{'group','session','day','nTrials','nCells','nGo','nNoGo','file'});
mf=fullfile(fileparts(OUTROOT),sprintf('manifest_%s.csv',GROUP));
writetable(M,mf);
fprintf('\n%d sessions in %.1f min\nmanifest -> %s\n', nS, toc(t0)/60, mf);
fprintf('DONE_MAKE_CACHE\n');


%% ================= helpers =================
function v=fastHex2double(hx)
    c=double(hx(3:end)); d=c-48; m=d>=49; d(m)=d(m)-39; m=d>=17; d(m)=d(m)-7;
    b=uint8(d(1:2:end)*16+d(2:2:end)); v=typecast(b(65:end),'double');
end

function row=writeSession(sessName,dayStr,trials,cells,vstim,traces,DBNAME,GROUP,OUTROOT, ...
        rasterLen,endIdx,bins,nWin,baselineBins,baselineWindow,SD_FLOOR,GLOB_FLOOR, ...
        centers,tvec,tvec_centre,dt,winFrames,stepFrames,tcWindowSec, ...
        cellMap,cellOSI,cellPref,trialMap,T2)
    row={};
    if isempty(trials), return; end
    uT=unique(trials); uC=unique(cells); nT=numel(uT); nC=numel(uC);
    if nT<5||nC<5, return; end

    Traw=NaN(nT,nC,rasterLen,'single');
    [~,lT]=ismember(trials,uT); [~,lC]=ismember(cells,uC);
    Vstim=NaN(nT,1);
    for r=1:numel(traces)
        x=fastHex2double(traces{r})'; x=x(1:min(end,rasterLen));
        if numel(x)<endIdx, continue; end
        Traw(lT(r),lC(r),1:numel(x))=single(x);
        Vstim(lT(r))=vstim(r);
    end

    % ---- raw window sums ----
    X=zeros(nT,nC,nWin);
    for w=1:nWin, X(:,:,w)=sum(Traw(:,:,bins{w}),3,'omitnan'); end
    assert(all(X(:)==round(X(:))) && max(X(:))<32767, 'window sums must be small integers');

    % ---- per-trial/cell baseline mean (for the subtraction variant) ----
    b=squeeze(mean(Traw(:,:,baselineWindow),3,'omitnan'));
    b(~isfinite(b))=0;

    % ---- normalisation constants ----
    Xb=reshape(permute(X(:,:,baselineBins),[1 3 2]),[],nC);       % baseline rows x cells
    mu_base=mean(Xb,1,'omitnan');
    sd_base_unfloored=std(Xb,0,1,'omitnan');
    sd_base=sd_base_unfloored;
    sd_base(~isfinite(sd_base))=0;
    sd_base=max(sd_base,SD_FLOOR);                                 % HARD floor 0.05

    Xsub=X-3*repmat(b,1,1,nWin);                                   % every bin is 3 frames
    Xs=reshape(permute(Xsub,[1 3 2]),[],nC);                       % all rows x cells
    mu_glob_sub=mean(Xs,1,'omitnan');
    sd_glob_sub=std(Xs,0,1,'omitnan');
    sd_glob_sub(sd_glob_sub<GLOB_FLOOR|~isfinite(sd_glob_sub))=GLOB_FLOOR;

    % ---- metadata ----
    cellID=uC(:); OSI=nan(nC,1); PrefOri=nan(nC,1);
    for i=1:nC
        key=[sessName '|' num2str(cellID(i))];
        if isKey(cellMap,key), j=cellMap(key); OSI(i)=cellOSI(j); PrefOri(i)=cellPref(j); end
    end
    trialID=uT(:); Response=nan(nT,1); Outcome=repmat({''},nT,1);
    for i=1:nT
        key=string([sessName '|' num2str(trialID(i))]);
        if isKey(trialMap,key), j=trialMap(key); Response(i)=T2.Response(j); Outcome{i}=T2.Outcome{j}; end
    end

    S=struct();
    S.X=int16(X); S.b=single(b);
    S.mu_base=mu_base; S.sd_base=sd_base; S.sd_base_unfloored=sd_base_unfloored;
    S.mu_glob_sub=mu_glob_sub; S.sd_glob_sub=sd_glob_sub;
    S.sd_floor=SD_FLOOR;
    S.trialID=trialID; S.Visual_Stim=Vstim; S.Response=Response; S.Outcome={Outcome{:}}';
    S.cellID=cellID; S.OSI=OSI; S.PrefOri=PrefOri;
    S.session=sessName; S.day=dayStr; S.group=GROUP; S.DB=DBNAME;
    S.bins=bins; S.centers=centers; S.tvec=tvec; S.tvec_centre=tvec_centre;
    S.baselineBins=baselineBins; S.baselineWindow=baselineWindow;
    S.dt=dt; S.rasterLen=rasterLen; S.winFrames=winFrames; S.stepFrames=stepFrames;
    S.tcWindowSec=tcWindowSec;
    S.provenance=struct('source','Zenodo-format CSV export (DatasetsCorbo2023)', ...
        'qc','PrefOri NOT NULL & isCell>0.8 & DeconvCorr>0.8 & spikeProb_TCTrials>0.1 & roundness>0.2 & nPix*magFactor>15', ...
        'blocks','Visual (Orientation Tuning excluded)', ...
        'windowing','backward, all bins full 3 frames: (c-2):c for c=9:3:84', ...
        'made',datestr(now,'yyyy-mm-dd HH:MM'));

    outDir=fullfile(OUTROOT,GROUP,dayStr);
    if ~exist(outDir,'dir'), mkdir(outDir); end
    f=fullfile(outDir,sprintf('%s_%s_win.mat',DBNAME,sessName));
    save(f,'-struct','S','-v7.3');
    row={GROUP,sessName,dayStr,nT,nC,sum(Vstim==45),sum(Vstim~=45 & isfinite(Vstim)),f};
end
