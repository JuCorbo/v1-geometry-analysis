%% fractionActive_pull.m  active-neuron fraction from the shared windowed cache
% Computes the binomial responsive-fraction from the shared intermediate via
% getWindowed('raw'). The 1 s windows are exact unions of 5 bins (baseline =
% bins 2-6, 0-1s = 7-11, 1-2s = 12-16, 2-3s = 17-21), so summing raw window
% sums over those bins reproduces the frame-level 1 s counts. A neuron counts
% as active when its window P(APrE>=1) exceeds the baseline Clopper-Pearson
% 95% upper CI. Writes cache/fractionActive_windowed.mat; the reference file
% fractionActive_pulled.mat (used by validate_windowed_cache.m) is left intact.
codeDir=fileparts(mfilename('fullpath')); aDir=fileparts(codeDir);
figDir=fileparts(fileparts(aDir)); sharedCode=fullfile(figDir,'_sharedWindowedData','code');
sharedData=fullfile(figDir,'_sharedWindowedData','data'); addpath(sharedCode);
cacheDir=fullfile(aDir,'cache'); if ~exist(cacheDir,'dir'), mkdir(cacheDir); end

behavCond={'D1','D2','D3','D4','D5','D6'};
noGoOri=[135,90,75,70,65,60]; goOri=45; gngAng=[90 45 30 25 20 15];
alpha_binom=0.05; nStimWins=3;
baseBins=2:6; winBins={7:11,12:16,17:21};                 % 1 s unions of 3-frame bins
condLab={'All','Go','NoGo'}; nCond=3;
colNames={'All_0_1','All_1_2','All_2_3','Go_0_1','Go_1_2','Go_2_3','NoGo_0_1','NoGo_1_2','NoGo_2_3'};
nCols=nCond*nStimWins;
grpFolders={'naive','trained'};
fracAct_binom=cell(2,numel(behavCond)); sessName_binom=cell(2,numel(behavCond));
for gi=1:2, for di=1:numel(behavCond), fracAct_binom{gi,di}=zeros(0,nCols); sessName_binom{gi,di}={}; end, end
NR=struct('group',[],'cond',[],'baseP',[],'baseUp',[],'N',[],'act',zeros(0,nStimWins), ...
          'pstim',zeros(0,nStimWins),'OSI',[],'dOri',[],'day',[]);

for gi=1:2
    for di=1:numel(behavCond)
        files=dir(fullfile(sharedData,grpFolders{gi},behavCond{di},'*_win.mat'));
        vGo=goOri; vNo=noGoOri(di);
        for si=1:numel(files)
            S=load(fullfile(files(si).folder,files(si).name));
            [Xr,info]=getWindowed(S,'raw');                 % [nTrials x nCells x nWin]
            V=info.Visual_Stim(:); osiC=info.OSI(:); prefC=info.PrefOri(:);
            nT=size(Xr,1); nC=size(Xr,2); if nT<5||nC<5, continue; end
            goM=(V==vGo); noM=(V==vNo); if ~any(goM)||~any(noM), continue; end
            baseCnt=squeeze(sum(Xr(:,:,baseBins),3,'omitnan'));       % nTrials x nCells (1 s baseline)
            stimCnt=cell(1,nStimWins);
            for w=1:nStimWins, stimCnt{w}=squeeze(sum(Xr(:,:,winBins{w}),3,'omitnan')); end
            masks={goM|noM, goM, noM}; row=NaN(1,nCols);
            for cc=1:nCond
                tm=masks{cc}; nTs=sum(tm);
                bC=baseCnt(tm,:); k=sum(bC>=1,1,'omitnan');
                ok=sum(isfinite(bC),1)>0; nOk=sum(ok); if nOk<3, continue; end
                ku=k(ok); up=betainv(1-alpha_binom/2,ku+1,nTs-ku); up(ku>=nTs)=1;
                actN=false(nOk,nStimWins); psN=nan(nOk,nStimWins);
                for w=1:nStimWins
                    sc=stimCnt{w}(tm,ok); ps=sum(sc>=1,1,'omitnan')/nTs;
                    row((cc-1)*nStimWins+w)=sum(ps>up)/nOk; actN(:,w)=(ps>up)'; psN(:,w)=ps';
                end
                if cc>=2
                    stimOri=(cc==2)*vGo+(cc==3)*vNo; dO=abs(prefC(ok)-stimOri); dO=min(dO,180-dO);
                    NR.group=[NR.group;repmat(gi,nOk,1)]; NR.cond=[NR.cond;repmat(cc,nOk,1)];
                    NR.baseP=[NR.baseP;(ku(:)/nTs)]; NR.baseUp=[NR.baseUp;up(:)]; NR.N=[NR.N;repmat(nTs,nOk,1)];
                    NR.act=[NR.act;actN]; NR.pstim=[NR.pstim;psN]; NR.OSI=[NR.OSI;osiC(ok)]; NR.dOri=[NR.dOri;dO(:)];
                    NR.day=[NR.day;repmat(di,nOk,1)];
                end
            end
            fracAct_binom{gi,di}=[fracAct_binom{gi,di};row];
            sessName_binom{gi,di}{end+1}=char(info.session);
        end
        fprintf('  %-7s | %s : %d sessions\n',grpFolders{gi},behavCond{di},size(fracAct_binom{gi,di},1));
    end
end
save(fullfile(cacheDir,'fractionActive_windowed.mat'),'fracAct_binom','sessName_binom', ...
     'behavCond','gngAng','colNames','condLab','NR','-v7.3');
fprintf('saved -> %s\n',fullfile(cacheDir,'fractionActive_windowed.mat'));
