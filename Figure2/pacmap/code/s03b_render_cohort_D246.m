%% Cohort PaCMAP render, even days (Fig S1)
% Same alignment/render pipeline as 03_render_cohort_and_theta_D135.m, plotting
% the even days D2/D4/D6: one 2x3 trajectory/cloud figure per DB, two-pass aligned.
% The theta angles, stats and example figure are produced by 03 (all days).

clear;

%% ---------------- CONFIG (GLOBAL) ----------------
baseRoot = getenv('PAC_BASEROOT');
if isempty(baseRoot)
    baseRoot = fullfile(fileparts(mfilename('fullpath')),'..','cache','main_pacmap_bwgz_pm080','baseline_seed0');
end

DBlist   = {'BatuhanADNaive', 'BatuhanAngDist'};   % naive, trained
DBpretty = containers.Map(DBlist, {'Naive','Trained'});

days       = 1:6;          % all days to process (alignment etc.)
daysToPlot = [2 4 6];      % subset to display, e.g. [1 3 5]

dimsToPlot   = [1 2 3];              % PaCMAP 3D
snapshotTime = 1.5;                  % cloud snapshot time (sec)
labelTimes      = [0, 0.5, 3];          % marker times for trajectories (sec)
showTimeMarkers = false;             % set true to show time markers on trajectories

% Angle computation (displacement vector method)
angleT0 = 0.0;      % seconds
angleT1 = 0.8;      % seconds
showAnglePoints = true;

% centroid smoothing (used when computing trajectory landmarks)
sgolayOrder = 2;
sgolayFrame = 3;                    % odd; auto-adjusted per session if needed

% Go/NoGo colors (within each DB supermouse fig)
colGo   = [0.20 0.40 0.80];
colNoGo = [0.80 0.20 0.20];

% DB comparison colors (combined angle fig)
colTrained = [0.20 0.40 0.80];      % blue
colNaive   = [0.85 0.65 0.13];      % gold-ish

% Styles
thinLW   = 0.8;
thickLW  = 3.0;
thinFade = 0.60;     % blend toward white for thin trajectories (color only)

ptSize   = 6;
padTraj  = 0.06;
padClouds= 0.30;

% Trajectory color timing behavior:
% gray before 0, colored from 0->3 with increasing intensity, gray after 3
trajGray = [0.60 0.60 0.60];
colorT0  = 0.0;
colorT1  = 3.0;
plotTmin = -1.0;     % clip trajectories before this time
plotTmax =  3.0;     % clip trajectories after this time

fadeGamma      = 1.01;  % >1 = stays faint longer, saturates later
fadeMinToWhite = 0.50;  % 0..1 : whiteness at start of colored epoch (higher = fainter at 0)

% Save
saveFigures = true;
savePdf = true;
saveFig = false;

% Output: supp panels, S1 names
outDirCombined = fullfile(fileparts(mfilename('fullpath')),'..','panels','supp');
nameCohort = 'S1_cohort_D246_';   % + Naive/Trained
if saveFigures && ~exist(outDirCombined, 'dir')
    mkdir(outDirCombined);
end

%% ---------------- RUN BOTH DBS ----------------
anglesAll  = struct();  % anglesAll.(DB).byDay{d} = [angles]
DayStruct  = struct();  % DayStruct.(DB) = Day array from processOneDB
Tall = table();        % long-format angle table built across DBs

for dbi = 1:numel(DBlist)
    DB = DBlist{dbi};
    rootDir = fullfile(baseRoot, DB);

    [anglesByDay, globalDayUsed, Tdb, DayStruct.(DB)] = processOneDB( ...
        rootDir, DB, DBpretty(DB), outDirCombined, nameCohort, ...
        days, daysToPlot, dimsToPlot, snapshotTime, labelTimes, ...
        angleT0, angleT1, ...
        sgolayOrder, sgolayFrame, ...
        colGo, colNoGo, ...
        thinLW, thickLW, thinFade, ...
        ptSize, padTraj, padClouds, ...
        trajGray, colorT0, colorT1, ...
        fadeGamma, fadeMinToWhite, ...
        plotTmin, plotTmax, ...
        showTimeMarkers, ...
        saveFigures, savePdf, saveFig, showAnglePoints);

    anglesAll.(DB) = struct('byDay',{anglesByDay}, 'globalDay', globalDayUsed);

    % append to global tall table
    Tall = [Tall; Tdb]; %#ok<AGROW>
end

%% =========================================================================
% LOCAL FUNCTION: process one DB (build Day struct, 2-pass align, plot 2x6, return anglesByDay + Tdb)
% =========================================================================
function [anglesByDay, globalDay, Tdb, DayOut] = processOneDB( ...
    rootDir, DBname, DBprettyName, outDir, nameCohort, ...
    days, daysToPlot, dimsToPlot, snapshotTime, labelTimes, ...
    angleT0, angleT1, ...
    sgolayOrder, sgolayFrame, ...
    colGo, colNoGo, ...
    thinLW, thickLW, thinFade, ...
    ptSize, padTraj, padClouds, ...
    trajGray, colorT0, colorT1, ...
    fadeGamma, fadeMinToWhite, ...
    plotTmin, plotTmax, ...
    showTimeMarkers, ...
    saveFigures, savePdf, saveFig, showAnglePoints)

    %#ok<INUSD> showAnglePoints

    % outDir and nameCohort are passed in
    if saveFigures && ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    % Long-format table for this DB
    angleDeg = [];
    dayCol   = [];
    training = strings(0,1);
    animalID = strings(0,1);
    sessPath = strings(0,1);

    % Gather files by day
    dayFiles = cell(1, max(days));
    for d = days
        dayFolder = fullfile(rootDir, sprintf('D%d', d));
        dayFiles{d} = dir(fullfile(dayFolder, '*_PaCMAP_UMATO.mat'));
    end

    % PASS 1 containers
    Day = repmat(struct( ...
        'exists', false, ...
        'refC', [], ...
        'tCentersWin', [], ...
        'sessions', [] ), 1, max(days));

    for d = days
        these = dayFiles{d};
        if isempty(these), continue; end

        SessRaw = struct([]);

        for i = 1:numel(these)
            matPath = fullfile(these(i).folder, these(i).name);
            S = load(matPath);

            if ~isfield(S,'Y_pacmap') || size(S.Y_pacmap,2) < max(dimsToPlot)
                continue;
            end
            req = {'labels','time_window_idx','time_window_center'};
            if any(~isfield(S, req))
                continue;
            end

            Y3     = S.Y_pacmap(:, dimsToPlot);
            labels = S.labels(:);
            twi    = S.time_window_idx(:);
            twc    = S.time_window_center(:);

            stimVals = unique(labels);
            if numel(stimVals) ~= 2
                continue;
            end

            % identify Go/NoGo (Go=45 if present)
            if any(stimVals == 45)
                goVal   = 45;
                nogoVal = stimVals(stimVals ~= 45); nogoVal = nogoVal(1);
            else
                goVal   = stimVals(1);
                nogoVal = stimVals(2);
            end

            uniqueTW = unique(twi);
            nWin     = numel(uniqueTW);

            tCentersWin = nan(nWin,1);
            for wi = 1:nWin
                wIdx = uniqueTW(wi);
                tCentersWin(wi) = mean(twc(twi==wIdx), 'omitnan');
            end

            % centroid trajectories
            trajGo = nan(nWin,3);
            trajNo = nan(nWin,3);
            for wi = 1:nWin
                wIdx = uniqueTW(wi);
                idxG = (labels==goVal)   & (twi==wIdx);
                idxN = (labels==nogoVal) & (twi==wIdx);
                if any(idxG), trajGo(wi,:) = mean(Y3(idxG,:),1,'omitnan'); end
                if any(idxN), trajNo(wi,:) = mean(Y3(idxN,:),1,'omitnan'); end
            end

            % smooth trajectories (optional; used for landmark stability)
            sgEff = sgolayFrame;
            if mod(sgEff,2)==0, sgEff = sgEff + 1; end
            if sgEff > nWin
                sgEff = max(3, 2*floor((nWin-1)/2)+1);
            end
            if sgEff >= 3 && sgEff <= nWin
                vg = all(isfinite(trajGo),2);
                vn = all(isfinite(trajNo),2);
                if sum(vg) >= sgEff
                    for dd2 = 1:3
                        trajGo(vg,dd2) = sgolayfilt(trajGo(vg,dd2), sgolayOrder, sgEff);
                    end
                end
                if sum(vn) >= sgEff
                    for dd2 = 1:3
                        trajNo(vn,dd2) = sgolayfilt(trajNo(vn,dd2), sgolayOrder, sgEff);
                    end
                end
            end

            Cstack = [trajGo; trajNo];  % landmarks (2*nWin x 3)

            % snapshot window index nearest to snapshotTime
            [~, wNearest] = min(abs(tCentersWin - snapshotTime));
            snapWIdx = uniqueTW(wNearest);
            snapT    = tCentersWin(wNearest);

            SessRaw(end+1).path = matPath; %#ok<AGROW>
            SessRaw(end).Y = Y3;
            SessRaw(end).labels = labels;
            SessRaw(end).twi = twi;
            SessRaw(end).goVal = goVal;
            SessRaw(end).nogoVal = nogoVal;
            SessRaw(end).uniqueTW = uniqueTW;
            SessRaw(end).tCentersWin = tCentersWin;
            SessRaw(end).trajGo = trajGo;
            SessRaw(end).trajNo = trajNo;
            SessRaw(end).Cstack = Cstack;
            SessRaw(end).snapWIdx = snapWIdx;
            SessRaw(end).snapT = snapT;
        end

        if isempty(SessRaw), continue; end

            % choose within-day reference session = medoid of Cstack landmarks
        nSess = numel(SessRaw);

        if nSess == 1
            refIdx = 1;
        else
            distMat = nan(nSess, nSess);

            for ii = 1:nSess
                Ci = SessRaw(ii).Cstack;

                for jj = 1:nSess
                    Cj = SessRaw(jj).Cstack;

                    good = all(isfinite(Ci),2) & all(isfinite(Cj),2);

                    if sum(good) < 10
                        distMat(ii,jj) = inf;
                    else
                        dvec = sqrt(sum((Ci(good,:) - Cj(good,:)).^2, 2));
                        distMat(ii,jj) = mean(dvec, 'omitnan');
                    end
                end
            end

            distScore = sum(distMat, 2, 'omitnan');

            % if some candidates are all inf/nan, fall back safely
            if all(~isfinite(distScore))
                refIdx = 1;
            else
                distScore(~isfinite(distScore)) = inf;
                [~, refIdx] = min(distScore);
            end
        end

        Cref_day = SessRaw(refIdx).Cstack;

        % align to within-day reference
        SessA = struct([]);
        for i = 1:numel(SessRaw)
            Cs   = SessRaw(i).Cstack;
            good = all(isfinite(Cref_day),2) & all(isfinite(Cs),2);

            if sum(good) < 10
                Ya      = SessRaw(i).Y;
                trajGoA = SessRaw(i).trajGo;
                trajNoA = SessRaw(i).trajNo;
            else
                [~, ~, tr] = procrustes(Cref_day(good,:), Cs(good,:), ...
                    'Scaling', false, 'Reflection','best');
                Trot = tr.T; b = tr.b; c = tr.c(1,:);

                Ya      = b * SessRaw(i).Y      * Trot + c;
                trajGoA = b * SessRaw(i).trajGo * Trot + c;
                trajNoA = b * SessRaw(i).trajNo * Trot + c;
            end

            idxSnap = (SessRaw(i).twi == SessRaw(i).snapWIdx);
            Ys      = Ya(idxSnap,:);
            lab     = SessRaw(i).labels(idxSnap);

            SessA(end+1).path = SessRaw(i).path; %#ok<AGROW>
            SessA(end).trajGo = trajGoA;
            SessA(end).trajNo = trajNoA;
            SessA(end).tCentersWin = SessRaw(i).tCentersWin;
            SessA(end).cloudY = Ys;
            SessA(end).cloudLabels = lab;
            SessA(end).goVal = SessRaw(i).goVal;
            SessA(end).nogoVal = SessRaw(i).nogoVal;
            SessA(end).snapT = SessRaw(i).snapT;
            SessA(end).fullY = Ya;                        % full aligned embedding (all windows)
            SessA(end).fullTwi = SessRaw(i).twi;         % window index per point
            SessA(end).fullLabels = SessRaw(i).labels;   % stimulus label per point
            SessA(end).uniqueTW = SessRaw(i).uniqueTW;   % ordered window indices
        end

        Day(d).exists      = true;
        Day(d).refC        = Cref_day;
        Day(d).tCentersWin = SessRaw(refIdx).tCentersWin;
        Day(d).sessions    = SessA;
    end

    % PASS 2: day -> global
    globalDay = 1;
    if ~(Day(globalDay).exists)
        globalDay = find([Day.exists], 1, 'first');
    end
    if isempty(globalDay)
        error('No usable days found under %s', rootDir);
    end

    Cref_global = Day(globalDay).refC;

    DayTr = repmat(struct('ok',false,'b',1,'T',eye(3),'c',[0 0 0]), 1, max(days));
    for d = days
        if ~Day(d).exists, continue; end
        Cref_day = Day(d).refC;

        good = all(isfinite(Cref_global),2) & all(isfinite(Cref_day),2);
        if sum(good) < 10
            warning('[%s] Too few landmarks to align day %d to global day %d; leaving as-is.', DBname, d, globalDay);
            continue;
        end

        [~, ~, trDay] = procrustes(Cref_global(good,:), Cref_day(good,:), ...
            'Scaling', false, 'Reflection','best');

        DayTr(d).ok = true;
        DayTr(d).b  = trDay.b;
        DayTr(d).T  = trDay.T;
        DayTr(d).c  = trDay.c(1,:);
    end

    for d = days
        if ~Day(d).exists || ~DayTr(d).ok, continue; end
        b    = DayTr(d).b;
        Trot = DayTr(d).T;
        c    = DayTr(d).c;

        for i = 1:numel(Day(d).sessions)
            Day(d).sessions(i).trajGo = b * Day(d).sessions(i).trajGo * Trot + c;
            Day(d).sessions(i).trajNo = b * Day(d).sessions(i).trajNo * Trot + c;
            Day(d).sessions(i).cloudY = b * Day(d).sessions(i).cloudY * Trot + c;
        end
    end

    % PLOT 2 x numel(daysToPlot) + compute angles
    nCols = numel(daysToPlot);
    figW  = max(600, 500 * nCols);
    fig = figure('Color','w','Position',[20 60 figW 950]);
    tl  = tiledlayout(2, nCols, 'TileSpacing','none','Padding','compact');
    sgtitle(tl, sprintf('%s PaCMAP | within-day aligned then day->D%d aligned', DBprettyName, globalDay), 'Interpreter','none');

    axTop = gobjects(1, nCols);
    axBot = gobjects(1, nCols);

    allTrajPts_GLOBAL  = [];
    allCloudPts_GLOBAL = [];

    colGoThin = blend_with_white(colGo, thinFade);
    colNoThin = blend_with_white(colNoGo, thinFade);

    anglesByDay = cell(1, max(days));

    % --- Pass A: compute angles for ALL days (needed for combined boxplot) ---
    for d = days
        if ~Day(d).exists, continue; end
        Sess = Day(d).sessions;
        for i = 1:numel(Sess)
            angDeg = angleBetweenDisplacements(Sess(i).trajGo, Sess(i).trajNo, Sess(i).tCentersWin, angleT0, angleT1);
            anglesByDay{d}(end+1,1) = angDeg; %#ok<AGROW>
            angleDeg(end+1,1) = angDeg;        %#ok<AGROW>
            dayCol(end+1,1)   = d;             %#ok<AGROW>
            training(end+1,1) = string(DBprettyName); %#ok<AGROW>
            sessPath(end+1,1) = string(Sess(i).path); %#ok<AGROW>
            animalID(end+1,1) = extractAnimalID_fromPath(Sess(i).path); %#ok<AGROW>
        end
    end

    % --- Pass B: plot only daysToPlot ---
    for di = 1:numel(daysToPlot)
        d = daysToPlot(di);
        axTop(di) = nexttile(tl, di);
        hold(axTop(di),'on'); grid(axTop(di),'on'); box(axTop(di),'on');
        title(axTop(di), sprintf('D%d trajectories', d), 'Interpreter','none');

        axBot(di) = nexttile(tl, nCols+di);
        hold(axBot(di),'on'); grid(axBot(di),'on'); box(axBot(di),'on');
        title(axBot(di), sprintf('D%d cloud @ %.2fs', d, snapshotTime), 'Interpreter','none');

        if ~Day(d).exists
            text(axTop(di), 0.5, 0.5, 'missing', 'Units','normalized', 'HorizontalAlignment','center');
            axis(axTop(di),'off');
            text(axBot(di), 0.5, 0.5, 'missing', 'Units','normalized', 'HorizontalAlignment','center');
            axis(axBot(di),'off');
            continue;
        end

        Sess = Day(d).sessions;
        if isempty(Sess)
            text(axTop(di), 0.5, 0.5, 'no usable sessions', 'Units','normalized', 'HorizontalAlignment','center');
            axis(axTop(di),'off');
            text(axBot(di), 0.5, 0.5, 'no usable sessions', 'Units','normalized', 'HorizontalAlignment','center');
            axis(axBot(di),'off');
            continue;
        end

        goTrajCell = cell(1,numel(Sess));
        noTrajCell = cell(1,numel(Sess));

        for i = 1:numel(Sess)
            tMaskThin = Sess(i).tCentersWin >= plotTmin & Sess(i).tCentersWin <= plotTmax;
            tgClip = Sess(i).trajGo; tgClip(~tMaskThin,:) = NaN;
            tnClip = Sess(i).trajNo; tnClip(~tMaskThin,:) = NaN;
            plot3(axTop(di), tgClip(:,1), tgClip(:,2), tgClip(:,3), '-', ...
                'Color', colGoThin, 'LineWidth', thinLW);
            plot3(axTop(di), tnClip(:,1), tnClip(:,2), tnClip(:,3), '-', ...
                'Color', colNoThin, 'LineWidth', thinLW);

            goTrajCell{i} = Sess(i).trajGo;
            noTrajCell{i} = Sess(i).trajNo;

            allTrajPts_GLOBAL = [allTrajPts_GLOBAL; Sess(i).trajGo; Sess(i).trajNo]; %#ok<AGROW>

            Ys  = Sess(i).cloudY;
            lab = Sess(i).cloudLabels;

            scatter3(axBot(di), Ys(lab==Sess(i).goVal,1),   Ys(lab==Sess(i).goVal,2),   Ys(lab==Sess(i).goVal,3), ...
                ptSize, colGo,   'filled','MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
            scatter3(axBot(di), Ys(lab==Sess(i).nogoVal,1), Ys(lab==Sess(i).nogoVal,2), Ys(lab==Sess(i).nogoVal,3),...
                ptSize, colNoGo, 'filled','MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);

            allCloudPts_GLOBAL = [allCloudPts_GLOBAL; Ys]; %#ok<AGROW>

            % (angle accumulation handled in Pass A above)
        end

        [meanGo, okG] = meanTraj(goTrajCell);
        [meanNo, okN] = meanTraj(noTrajCell);
        if ~(okG && okN), continue; end

        tWin = Sess(1).tCentersWin;
        plotTimeGrayColorGrayTrajectory(axTop(di), meanGo, tWin, colGo,   thickLW, colorT0, colorT1, trajGray, fadeGamma, fadeMinToWhite, plotTmin, plotTmax);
        plotTimeGrayColorGrayTrajectory(axTop(di), meanNo, tWin, colNoGo, thickLW, colorT0, colorT1, trajGray, fadeGamma, fadeMinToWhite, plotTmin, plotTmax);

        markerList  = {'o','^','s'};
        if showTimeMarkers
            labelTimesClip = labelTimes(labelTimes >= plotTmin & labelTimes <= plotTmax);
            addTimeMarkersOnTrajectory(axTop(di), meanGo, tWin, labelTimesClip, markerList, colGo);
            addTimeMarkersOnTrajectory(axTop(di), meanNo, tWin, labelTimesClip, markerList, colNoGo);
        end

        xlabel(axTop(di),'dim1'); ylabel(axTop(di),'dim2'); zlabel(axTop(di),'dim3');
        xlabel(axBot(di),'dim1'); ylabel(axBot(di),'dim2'); zlabel(axBot(di),'dim3');
        view(axTop(di),3); view(axBot(di),3);
    end

    trajLims_GLOBAL  = computeLims(allTrajPts_GLOBAL,  padTraj);
    cloudLims_GLOBAL = computeLims(allCloudPts_GLOBAL, padClouds);

    % --- Compute optimal camera angle from D1 (globalDay) aligned trajectories ---
    % Use the plane spanned by:
    %   v1 = mean Go displacement from angleT0 to angleT1  (trajectory direction)
    %   v2 = Go-NoGo centroid separation at angleT1        (divergence direction)
    % Camera looks along the normal to this plane; up = v1 (trajectory direction).
    camPos = [];   % fallback: empty = use MATLAB default view(3)
    camUp  = [];
    refDay = globalDay;
    if Day(refDay).exists && ~isempty(Day(refDay).sessions)
        SessRef = Day(refDay).sessions;
        goCell = {SessRef.trajGo};
        noCell = {SessRef.trajNo};
        [mGref, okG] = meanTraj(goCell);
        [mNref, okN] = meanTraj(noCell);
        tRef = SessRef(1).tCentersWin(:);

        if okG && okN
            [tRefS, ord] = sort(tRef);
            mGref = mGref(ord,:);
            mNref = mNref(ord,:);

            % Interpolate to angleT0 and angleT1
            p0G = interp1(tRefS, mGref, angleT0, 'pchip', 'extrap');
            p1G = interp1(tRefS, mGref, angleT1, 'pchip', 'extrap');
            p0N = interp1(tRefS, mNref, angleT0, 'pchip', 'extrap');
            p1N = interp1(tRefS, mNref, angleT1, 'pchip', 'extrap');

            v1 = p1G - p0G;                  % Go trajectory direction
            v2 = (p1G + p1N)/2 - (p0G + p0N)/2; % mean displacement (both stims)
            vSep = p1G - p1N;                % Go-NoGo separation at angleT1

            % Use v1 (traj direction) as up, vSep as the other spanning vector
            if norm(v1) > eps && norm(vSep) > eps
                v1   = v1   / norm(v1);
                vSep = vSep / norm(vSep);
                % Remove component of vSep along v1
                vSep = vSep - dot(vSep, v1) * v1;
                if norm(vSep) > eps
                    vSep   = vSep / norm(vSep);
                    normal = cross(v1, vSep);   % look-from direction
                    normal = normal / norm(normal);

                    % Camera sits along the normal, looking at the scene center
                    sceneCtr = mean(trajLims_GLOBAL, 2)';        % [1x3]
                    camDist  = 3 * max(trajLims_GLOBAL(:,2) - trajLims_GLOBAL(:,1));
                    normal   = normal(:)';                        % enforce [1x3]
                    camPos   = sceneCtr + camDist * normal;       % [1x3]
                    camUp    = v1(:)';                            % enforce [1x3]
                end
            end
        end
    end

    for di = 1:numel(daysToPlot)
        if isgraphics(axTop(di))
            setAxes3(axTop(di), trajLims_GLOBAL);
            grid(axTop(di),'on'); box(axTop(di),'on');
            if ~isempty(camPos)
                set(axTop(di), 'CameraPosition', camPos, 'CameraUpVector', camUp, 'CameraViewAngleMode','auto');
            else
                view(axTop(di),3);
            end
        end
        if isgraphics(axBot(di))
            setAxes3(axBot(di), cloudLims_GLOBAL);
            grid(axBot(di),'on'); box(axBot(di),'on');
            if ~isempty(camPos)
                set(axBot(di), 'CameraPosition', camPos, 'CameraUpVector', camUp, 'CameraViewAngleMode','auto');
            else
                view(axBot(di),3);
            end
        end
    end

    rotate3d(fig,'on');
    hTop = linkprop(axTop(isgraphics(axTop)), {'CameraPosition','CameraTarget','CameraUpVector','CameraViewAngle'});
    hBot = linkprop(axBot(isgraphics(axBot)), {'CameraPosition','CameraTarget','CameraUpVector','CameraUpVector','CameraViewAngle'});
    setappdata(fig, 'CamLinkTop', hTop);
    setappdata(fig, 'CamLinkBot', hBot);

    if saveFigures
        outBase = [nameCohort char(DBprettyName)];
        if savePdf
            exportgraphics(fig, fullfile(outDir, [outBase '.pdf']), 'ContentType','vector', 'BackgroundColor','white');
            exportgraphics(fig, fullfile(outDir, [outBase '.png']), 'Resolution',140, 'BackgroundColor','white');
        end
        if saveFig
            savefig(fig, fullfile(outDir, [outBase '.fig']));
        end
        fprintf('[%s] Saved main figure to %s\n', DBname, outDir);
    end

    % Export Day struct for downstream example figure
    DayOut = Day;

    % Build output table for this DB
    Tdb = table(angleDeg, dayCol, training, animalID, sessPath, ...
        'VariableNames', {'angleDeg','day','training','animalID','sessPath'});
end

function lims = computeLims(P, pad)
    if isempty(P)
        lims = [0 1; 0 1; 0 1];
        return;
    end
    mins = min(P,[],1);
    maxs = max(P,[],1);
    rngs = maxs - mins;
    mins = mins - (pad*rngs);
    maxs = maxs + (pad*rngs);
    lims = [mins(1) maxs(1); mins(2) maxs(2); mins(3) maxs(3)];
end

function setAxes3(ax, lims)
    xlim(ax, lims(1,:));
    ylim(ax, lims(2,:));
    zlim(ax, lims(3,:));
    axis(ax,'vis3d');
end

function c = blend_with_white(baseCol, fracWhite)
    fracWhite = min(max(fracWhite,0),1);
    c = (1-fracWhite)*baseCol + fracWhite*[1 1 1];
end

function [m, ok] = meanTraj(trajCell)
    ok = false; m = [];
    if isempty(trajCell), return; end
    nWin = size(trajCell{1},1);
    M = nan(nWin,3,numel(trajCell));
    for i = 1:numel(trajCell)
        T = trajCell{i};
        if size(T,1) ~= nWin || size(T,2) ~= 3
            return;
        end
        M(:,:,i) = T;
    end
    m = mean(M, 3, 'omitnan');
    ok = true;
end

function plotTimeGrayColorGrayTrajectory(ax, trajXYZ, tCentersWin, baseCol, lw, t0, t1, grayCol, gamma, minToWhite, plotTmin, plotTmax)
    if nargin < 12, plotTmin = -inf; plotTmax = inf; end
    if isempty(trajXYZ) || size(trajXYZ,1) < 2, return; end
    if numel(tCentersWin) ~= size(trajXYZ,1)
        tCentersWin = (1:size(trajXYZ,1))';
    end
    % Clip to plot window
    tMask = tCentersWin >= plotTmin & tCentersWin <= plotTmax;
    trajXYZ(~tMask,:) = NaN;
    t = tCentersWin(:);
    nWin = size(trajXYZ,1);

    % Flat color: draw entire clipped trajectory in baseCol
    % (per-segment gradient logic preserved below for reference)
    plot3(ax, trajXYZ(:,1), trajXYZ(:,2), trajXYZ(:,3), '-', 'Color', baseCol, 'LineWidth', lw);

    % --- gradient logic (kept for reference, not executed) ---
    % for w = 1:(nWin-1)
    %     seg = trajXYZ(w:w+1,:);
    %     if ~all(isfinite(seg), 'all'), continue; end
    %     tMid = mean(t(w:w+1), 'omitnan');
    %     if ~isfinite(tMid) || tMid < t0 || tMid > t1
    %         segCol = grayCol;
    %     else
    %         frac = (tMid - t0) / max(eps, (t1 - t0));
    %         frac = min(max(frac, 0), 1);
    %         frac = frac.^gamma;
    %         whiteFrac = minToWhite * (1 - frac);
    %         segCol    = whiteFrac*[1 1 1] + (1-whiteFrac)*baseCol;
    %     end
    %     plot3(ax, seg(:,1), seg(:,2), seg(:,3), '-', 'Color', segCol, 'LineWidth', lw);
    % end
end

function addTimeMarkersOnTrajectory(ax, trajXYZ, tCentersWin, timesToMark, markerList, baseCol)
    if isempty(trajXYZ) || size(trajXYZ,1) < 2 || numel(tCentersWin) ~= size(trajXYZ,1)
        return;
    end

    t = tCentersWin(:);
    X = trajXYZ(:,1); Y = trajXYZ(:,2); Z = trajXYZ(:,3);
    [t, ord] = sort(t);
    X = X(ord); Y = Y(ord); Z = Z(ord);

    [tUnique, ~, ic] = unique(t);
    if numel(tUnique) < numel(t)
        X = accumarray(ic, X, [], @mean);
        Y = accumarray(ic, Y, [], @mean);
        Z = accumarray(ic, Z, [], @mean);
        t = tUnique;
    end

    timesToMark = timesToMark(:);
    xq = interp1(t, X, timesToMark, 'pchip', 'extrap');
    yq = interp1(t, Y, timesToMark, 'pchip', 'extrap');
    zq = interp1(t, Z, timesToMark, 'pchip', 'extrap');

    markerSize = 5;
    edgeLW     = 1.0;

    for k = 1:numel(timesToMark)
        mk = markerList{min(k, numel(markerList))};
        scatter3(ax, xq(k), yq(k), zq(k), markerSize^2, ...
            'Marker', mk, ...
            'MarkerFaceColor', baseCol, ...
            'MarkerEdgeColor', baseCol, ...
            'LineWidth', edgeLW);
    end
end

function angDeg = angleBetweenDisplacements(trajGo, trajNo, tCentersWin, t0, t1)
    if isempty(trajGo) || isempty(trajNo) || size(trajGo,2) ~= 3 || size(trajNo,2) ~= 3
        angDeg = nan; return;
    end
    if numel(tCentersWin) ~= size(trajGo,1) || numel(tCentersWin) ~= size(trajNo,1)
        angDeg = nan; return;
    end

    t = tCentersWin(:);
    [tS, ord] = sort(t);
    G = trajGo(ord,:);
    N = trajNo(ord,:);

    [tU, ~, ic] = unique(tS);
    if numel(tU) < numel(tS)
        G = [accumarray(ic, G(:,1), [], @mean), accumarray(ic, G(:,2), [], @mean), accumarray(ic, G(:,3), [], @mean)];
        N = [accumarray(ic, N(:,1), [], @mean), accumarray(ic, N(:,2), [], @mean), accumarray(ic, N(:,3), [], @mean)];
        tS = tU;
    end

    p0G = [interp1(tS, G(:,1), t0, 'pchip', nan), interp1(tS, G(:,2), t0, 'pchip', nan), interp1(tS, G(:,3), t0, 'pchip', nan)];
    p1G = [interp1(tS, G(:,1), t1, 'pchip', nan), interp1(tS, G(:,2), t1, 'pchip', nan), interp1(tS, G(:,3), t1, 'pchip', nan)];
    p0N = [interp1(tS, N(:,1), t0, 'pchip', nan), interp1(tS, N(:,2), t0, 'pchip', nan), interp1(tS, N(:,3), t0, 'pchip', nan)];
    p1N = [interp1(tS, N(:,1), t1, 'pchip', nan), interp1(tS, N(:,2), t1, 'pchip', nan), interp1(tS, N(:,3), t1, 'pchip', nan)];

    if any(~isfinite([p0G p1G p0N p1N]))
        angDeg = nan; return;
    end

    vG = p1G - p0G;
    vN = p1N - p0N;

    if norm(vG) < eps || norm(vN) < eps
        angDeg = nan; return;
    end

    cosang = dot(vG, vN) / (norm(vG) * norm(vN));
    cosang = min(max(cosang, -1), 1);
    angDeg = acosd(cosang);
end

function id = extractAnimalID_fromPath(pth)
% Conservative animal ID parsing from a session file path.
% Looks for tokens like 'B10_111020' or 'C02_121120', i.e. letters+digits '_' digits.
% If not found, uses basename minus common suffixes.

    s = char(pth);

    tok = regexp(s, '([A-Za-z]+\d+_\d+)', 'tokens', 'once');
    if ~isempty(tok)
        id = string(tok{1});
        return;
    end

    [~, b, ~] = fileparts(s);
    b = regexprep(b, '_PaCMAP_UMATO$', '');
    b = regexprep(b, '_forDR$', '');
    id = string(b);
end

