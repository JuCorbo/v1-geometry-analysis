%% runPupilCapacityFull.m
% Launcher for the full pupil capacity run (resumable).
addpath(fileparts(mfilename('fullpath')));
RUN_SESSION_MAX    = Inf;
RUN_CENTERS_SUBSET = [];
RUN_WINDOWS        = {'baseline','post'};
RUN_RESUME         = true;   % skip partitions already computed
fprintf('PUPIL_FULLRUN_START %s\n', datestr(now));
tAll = tic;
pupilCapacitySplits;
fprintf('PUPIL_FULLRUN_DONE %s  ELAPSED_HR=%.2f\n', datestr(now), toc(tAll)/3600);
