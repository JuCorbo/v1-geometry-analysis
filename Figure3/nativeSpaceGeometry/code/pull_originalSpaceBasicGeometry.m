%% pull_originalSpaceBasicGeometry.m
% Native-space Go/NoGo geometry per session and time window, RMS-normalized (/sqrt(nCells)):
%   centDist  : Euclidean distance between Go and NoGo centroids
%   radiusMix : mean single-trial distance to class centroid, averaged over Go/NoGo
%   dispMix   : std of single-trial distances to centroid, averaged over Go/NoGo
%   radiusGo/No, dispGo/No, centMagGo/No : class-specific
% Reads getWindowed('baseline') (baseline z-score, sigma floor 0.12). No SQL.
% Saves originalSpaceBasicGeometrySession.mat with the cell structure
% {group (1=naive, 2=trained), day (1..6)}{session} that make_originalSpaceGeometryPlotting.m consumes.
clear; clc;
codeDir    = fileparts(mfilename('fullpath'));
aDir       = fileparts(codeDir);                       % nativeSpaceGeometry
figDir     = fileparts(fileparts(aDir));               % Figures_final
sharedCode = fullfile(figDir,'_sharedWindowedData','code');
sharedData = fullfile(figDir,'_sharedWindowedData','data');
addpath(sharedCode);

DBlist    = {'BatuhanADNaive','BatuhanAngDist'};       % naive(1), trained(2)
groupDir  = {'naive','trained'};
BehavCond = {'D1','D2','D3','D4','D5','D6'};
GoOri     = 45;
% style vars (consumed by make)
colNaive=[0.85 0.65 0.13]; colTrain=[0.20 0.40 0.80]; thinAlpha=0.25; thinLW=0.8; thickLW=3.0; %#ok<NASGU>

nDay = numel(BehavCond);
[centDist_bySess,radius_bySess,disp_bySess,radiusGo_bySess,radiusNo_bySess, ...
 dispGo_bySess,dispNo_bySess,centMagGo_bySess,centMagNo_bySess,sessNames_byDay] = deal(cell(2,nDay));
tvec = [];

for gi = 1:2
    for d = 1:nDay
        files = dir(fullfile(sharedData, groupDir{gi}, sprintf('D%d',d), '*_win.mat'));
        for k = 1:numel(files)
            S = load(fullfile(files(k).folder, files(k).name));
            [X, info] = getWindowed(S, 'baseline');    % [nTr x nC x nW]
            if isempty(tvec), tvec = info.tvec(:)'; end
            vs = double(info.Visual_Stim(:));
            goM = (vs == GoOri); noM = (vs ~= GoOri & isfinite(vs));   % NoGo = day-specific other stimulus
            if sum(goM) < 2 || sum(noM) < 2, continue; end
            nW = size(X,3); nC = size(X,2); rmsDen = sqrt(nC);
            [centDist,radiusMix,dispMix,radiusGo,radiusNo,dispGo,dispNo,centMagGo,centMagNo] = deal(nan(1,nW));
            for w = 1:nW
                Xw = X(:,:,w); Xg = Xw(goM,:); Xn = Xw(noM,:);
                Xg = Xg(all(isfinite(Xg),2),:); Xn = Xn(all(isfinite(Xn),2),:);
                if size(Xg,1) < 2 || size(Xn,1) < 2, continue; end
                cg = mean(Xg,1); cn = mean(Xn,1);
                centDist(w)  = norm(cg - cn) / rmsDen;
                centMagGo(w) = norm(cg) / rmsDen; centMagNo(w) = norm(cn) / rmsDen;
                dg = vecnorm(Xg - cg, 2, 2) / rmsDen; dn = vecnorm(Xn - cn, 2, 2) / rmsDen;
                radiusGo(w) = mean(dg,'omitnan'); radiusNo(w) = mean(dn,'omitnan');
                dispGo(w)   = std(dg,0,'omitnan'); dispNo(w)  = std(dn,0,'omitnan');
                radiusMix(w) = mean([radiusGo(w) radiusNo(w)],'omitnan');
                dispMix(w)   = mean([dispGo(w) dispNo(w)],'omitnan');
            end
            centDist_bySess{gi,d}{end+1} = centDist;  radius_bySess{gi,d}{end+1} = radiusMix;  disp_bySess{gi,d}{end+1} = dispMix; %#ok<*SAGROW>
            radiusGo_bySess{gi,d}{end+1} = radiusGo;  radiusNo_bySess{gi,d}{end+1} = radiusNo;
            dispGo_bySess{gi,d}{end+1}   = dispGo;    dispNo_bySess{gi,d}{end+1}   = dispNo;
            centMagGo_bySess{gi,d}{end+1}= centMagGo; centMagNo_bySess{gi,d}{end+1}= centMagNo;
            sessNames_byDay{gi,d}{end+1} = char(info.session);
        end
        fprintf('  %-7s D%d: %d sessions\n', groupDir{gi}, d, numel(sessNames_byDay{gi,d}));
    end
end

outMat = fullfile(aDir,'cache','originalSpaceBasicGeometrySession.mat');
if ~exist(fullfile(aDir,'cache'),'dir'), mkdir(fullfile(aDir,'cache')); end
save(outMat, 'centDist_bySess','radius_bySess','disp_bySess', ...
    'radiusGo_bySess','radiusNo_bySess','dispGo_bySess','dispNo_bySess', ...
    'centMagGo_bySess','centMagNo_bySess','sessNames_byDay','tvec','BehavCond','DBlist', ...
    'colNaive','colTrain','thinAlpha','thinLW','thickLW','-v7.3');
fprintf('saved %s\n', outMat);
