%% 01b_pull_data4timeDR_baselinez.m
% SUPP Fig 2 (S2) + multiDR (S3) input: baseline z-score (sigma floor 0.12,
% applied by getWindowed) with NO per-trial subtraction, via getWindowed('baseline').
% Writes <DB>_<session>_TimeDR_forDR.mat; X_raw keeps the raw window sums ('raw').
codeDir = fileparts(mfilename('fullpath'));
outRoot = getenv('F2_BASELINE_OUT');
if isempty(outRoot)
    outRoot = fullfile(tempdir,'f2_forDR','baseline');
end
pull_forDR_fromCache('baseline', outRoot);
