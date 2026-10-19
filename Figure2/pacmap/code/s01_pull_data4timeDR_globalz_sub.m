%% 01_pull_data4timeDR_globalz_sub.m
% MAIN Fig 2 input: per-trial baseline subtracted then GLOBAL per-neuron z-score,
% via getWindowed('global'). Writes <DB>_<session>_TimeDR_forDR.mat for the PaCMAP/
% UMATO embedder. X_raw keeps the pre-z baseline-subtracted window sums ('sub').
% Output goes to a neutral master dir; the run step embeds it per seed.
codeDir = fileparts(mfilename('fullpath'));
outRoot = getenv('F2_GLOBAL_OUT');
if isempty(outRoot)
    outRoot = fullfile(tempdir,'f2_forDR','global');
end
pull_forDR_fromCache('global', outRoot);
