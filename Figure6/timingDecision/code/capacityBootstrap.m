%% Bootstrap CI for capacity divergence onset; merge into onsetCI.mat
clear; clc;
codeDir = fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
anaDir = fileparts(codeDir);
outDir = codeDir;
addpath(fullfile(fileparts(fileparts(anaDir)),'_sharedWindowedData','code'));  % groupSeparationTest
rng(7);
tcDir = fullfile(anaDir,'cache','capacityTimecourse');

loadGrp = @(db) local_load(tcDir, db);
[N,tvec] = loadGrp('BatuhanADNaive');    % [42 x 26]
[Tt,~]   = loadGrp('BatuhanAngDist');    % [60 x 26]
fprintf('loaded capacity: naive %dx%d, trained %dx%d\n', size(N), size(Tt));

% capacity params (match rollingCapacityPlottingPlusTiming): tail right, minSustained 2
onsetOffsetSec = 1*0.0646;  % 1*dt: report at window centre (tvec_centre)
gargs = {'tail','right','correction','cluster','alpha',0.05,'minSustained',2, ...
         'minHedgesG',0.2,'tStart',0,'onsetOffsetSec',onsetOffsetSec,'verbose',false};

rPt = groupSeparationTest(N,Tt,tvec,gargs{:},'nPerms',10000,'seed',1);
fprintf('point onset = %.3f s\n', rPt.onsetTime);

Nboot=1000; nPerms=1000; nN=size(N,1); nT=size(Tt,1);
onB=nan(Nboot,1);
for b=1:Nboot
    iN=randi(nN,nN,1); iT=randi(nT,nT,1);
    rb=groupSeparationTest(N(iN,:),Tt(iT,:),tvec,gargs{:},'nPerms',nPerms,'seed',b);
    onB(b)=rb.onsetTime;
end
onBd=onB(isfinite(onB)); ci=prctile(onBd,[2.5 97.5]); detRate=mean(isfinite(onB));
fprintf('capacity  | onset=%.2fs | 95%% CI [%.2f, %.2f] | bootMed=%.2f | detected %.0f%% | nN=%d nT=%d\n', ...
    rPt.onsetTime, ci(1), ci(2), median(onBd), 100*detRate, nN, nT);

% merge into onsetCI.mat
R=load(fullfile(anaDir,'data','onsetCI.mat')); results=R.results;
results.capacity = struct('onsetPoint',rPt.onsetTime,'ci',ci,'median',median(onBd), ...
    'detRate',detRate,'boot',onB,'nNaive',nN,'nTrained',nT,'tail','right');
save(fullfile(anaDir,'data','onsetCI.mat'),'results');
fprintf('merged capacity into onsetCI.mat\n');

function [M,tvec]=local_load(d,db)
    f=dir(fullfile(d,sprintf('*_%s_capacity_timecourse.mat',db)));
    f=f(~contains({f.name},'randThirdAvg'));
    M=[]; tvec=[];
    for i=1:numel(f)
        S=load(fullfile(d,f(i).name),'capacity_rolling','tvec');
        if isempty(tvec), tvec=S.tvec(:)'; end
        M(end+1,:)=S.capacity_rolling(:)'; %#ok<AGROW>
    end
end
