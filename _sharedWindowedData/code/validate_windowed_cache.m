%% validate_windowed_cache.m -- prove the intermediate reproduces SQL ground truth
%
% Recomputes Figure 1's fraction-of-active-neurons from the windowed cache and
% compares it, session by session, to the values pulled from the live SQL DB
% (cache/fractionActive_pulled.mat).  That metric uses 1 s windows which, under the
% adopted binning, are exact unions of 5 consecutive bins:
%     baseline frames 10-24 = bins  2-6
%     0-1 s    frames 25-39 = bins  7-11
%     1-2 s    frames 40-54 = bins 12-16
%     2-3 s    frames 55-69 = bins 17-21
% so the comparison should be EXACT. It validates decode, frame alignment, QC
% inclusion, trial labelling and binning in one shot -- including for naive, which
% has no other independent reference.

if ~exist('GROUP','var'), GROUP='trained'; end
codeDir=fileparts(mfilename('fullpath')); addpath(codeDir);          % getWindowed
DATAROOT=fullfile(fileparts(codeDir),'data');
FIG1=fullfile(fileparts(fileparts(codeDir)),'Figure1','sparsening');

C=load(fullfile(FIG1,'cache','fractionActive_pulled.mat'));
gi = 1 + strcmp(GROUP,'trained');            % 1 = naive, 2 = trained
behavCond=C.behavCond; noGoOri=[135,90,75,70,65,60]; goOri=45;
alpha_binom=0.05;
baseBins=2:6; stimBinSets={7:11, 12:16, 17:21};

d=dir(fullfile(DATAROOT,GROUP,'**','*_win.mat'));
fprintf('%s: %d session files\n', GROUP, numel(d));

got=containers.Map();
for i=1:numel(d)
    S=load(fullfile(d(i).folder,d(i).name));
    X=getWindowed(S,'raw');                       % raw integer sums
    di=find(strcmp(behavCond,S.day));
    V=S.Visual_Stim(:);
    goM=(V==goOri); noM=(V==noGoOri(di));
    if ~any(goM)||~any(noM), continue; end
    masks={goM|noM, goM, noM};
    row=NaN(1,9);
    for cc=1:3
        tm=masks{cc}; nTs=sum(tm);
        bC=sum(X(tm,:,baseBins),3);               % frames 10..24
        k=sum(bC>=1,1);
        nOk=size(X,2); if nOk<3, continue; end
        up=betainv(1-alpha_binom/2,k+1,nTs-k); up(k>=nTs)=1;
        for w=1:3
            sc=sum(X(tm,:,stimBinSets{w}),3);
            ps=sum(sc>=1,1)/nTs;
            row((cc-1)*3+w)=sum(ps>up)/nOk;
        end
    end
    got(S.session)=row;
end

fprintf('\n%-14s %-4s | %-26s | %-26s | maxdiff\n','session','day','cache (All 0-1 1-2 2-3)','SQL (All 0-1 1-2 2-3)');
fprintf('%s\n',repmat('-',1,96));
alld=[]; nMatch=0; nComp=0; missing={};
for di=1:numel(behavCond)
    names=C.sessName_binom{gi,di}; vals=C.fracAct_binom{gi,di};
    for k2=1:numel(names)
        sn=names{k2};
        if ~isKey(got,sn), missing{end+1}=sn; continue; end %#ok<AGROW>
        a=got(sn); b=vals(k2,:); dd=max(abs(a-b));
        alld(end+1)=dd; nComp=nComp+1; if dd<1e-12, nMatch=nMatch+1; end %#ok<AGROW>
        fprintf('%-14s %-4s | %7.4f %7.4f %7.4f | %7.4f %7.4f %7.4f | %.2e\n', ...
            sn, behavCond{di}, a(1),a(2),a(3), b(1),b(2),b(3), dd);
    end
end
fprintf('\n%d sessions compared; %d EXACT; max diff = %.3e\n', nComp, nMatch, max([alld 0]));
if ~isempty(missing), fprintf('missing from cache: %s\n', strjoin(missing,', ')); end
if nComp>0 && nMatch==nComp
    fprintf('\n*** %s VALIDATED against the live-DB values ***\n', upper(GROUP));
end

%% population mean under the adopted conventions, for reference
G=[];N=[];
for i=1:numel(d)
    S=load(fullfile(d(i).folder,d(i).name));
    Z=getWindowed(S,'baseline');                  % baseline z, sigma floor 0.12
    di=find(strcmp(behavCond,S.day)); V=S.Visual_Stim(:);
    g=V==goOri; n=V==noGoOri(di);
    if ~any(g)||~any(n), continue; end
    G(end+1,:)=squeeze(mean(mean(Z(g,:,:),1),2))'; %#ok<AGROW>
    N(end+1,:)=squeeze(mean(mean(Z(n,:,:),1),2))'; %#ok<AGROW>
end
S1=load(fullfile(d(1).folder,d(1).name),'tvec','tvec_centre');
fprintf('\npop-mean (%s, n=%d, backward full bins, sigma floor 0.12)\n', GROUP, size(G,1));
fprintf('  baseline bins 1-6 Go mean = %+.5f  (should be ~0)\n', mean(mean(G(:,1:6))));
[pk,ip]=max(mean(G,1));
fprintf('  Go peak %.4f at t_end=%.4f s (centre %.4f s)\n', pk, S1.tvec(ip), S1.tvec_centre(ip));
fprintf('\nDONE_VALIDATE\n');
