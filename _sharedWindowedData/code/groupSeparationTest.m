function results = groupSeparationTest(dataA, dataB, tvec, varargin)
% groupSeparationTest  Cluster-based permutation test (Maris & Oostenveld 2007),
% unpaired (between-group, label shuffle) or paired (within-group, sign-flip).
%
% Shared copy used by all figures (addpath _sharedWindowedData/code), so the
% cluster logic is defined in exactly one place.
%
% USAGE
%   results = groupSeparationTest(dataA, dataB, tvec, Name, Value, ...)
%
% INPUTS
%   dataA, dataB   [nA x nTime], [nB x nTime]  one row per session
%   tvec           [1 x nTime] time axis (seconds)
%   meanDiff is mean(dataB) - mean(dataA) (positive = B > A).
%
% KEY OPTIONS
%   'paired'      false (default). If true, nA must equal nB and rows must be
%                 aligned (same sessions); the null is a per-session SIGN-FLIP of
%                 the paired difference d = dataB - dataA, and the effect size is
%                 the paired Hedges g (standardised by SD of d). Use for
%                 within-group Go-vs-NoGo. Pass (Go, NoGo) => meanDiff = NoGo - Go.
%   'correction'  'cluster' (default) or 'bh'
%   'nPerms'      10000     'alpha' 0.05     'primaryAlpha' 0.05
%   'tail'        'both' (default), 'right', 'left'
%   'minSustained' 2        'minHedgesG' 0.1     'onsetOffsetSec' 0     'tStart' 0
%   'measureName' ''        'seed' 42            'verbose' true
%
% OUTPUT fields: pvals, qvals, clusterMask, clusterMass, clusterP, clusterPmin,
%   sigMask, sigMaskSust, onsetTime, hedgesG, meanDiff, meanA, meanB, tvec,
%   measureName, params.

%% ---- Input parsing --------------------------------------------------------
p = inputParser;
addRequired(p, 'dataA',   @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'dataB',   @(x) isnumeric(x) && ismatrix(x));
addRequired(p, 'tvec',    @(x) isnumeric(x) && isvector(x));
addParameter(p, 'paired',        false,  @islogical);
addParameter(p, 'correction',   'cluster', @(x) ismember(x, {'cluster','bh'}));
addParameter(p, 'nPerms',       10000,  @(x) isnumeric(x) && x > 0);
addParameter(p, 'alpha',        0.05,   @(x) isnumeric(x) && x > 0 && x < 1);
addParameter(p, 'primaryAlpha', 0.05,   @(x) isnumeric(x) && x > 0 && x < 1);
addParameter(p, 'tail',         'both', @(x) ismember(x, {'both','left','right'}));
addParameter(p, 'minSustained', 2,      @(x) isnumeric(x) && x >= 1);
addParameter(p, 'onsetOffsetSec', 0,    @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'tStart',       0,      @isnumeric);
addParameter(p, 'minHedgesG',   0.1,    @(x) isnumeric(x) && x >= 0);
addParameter(p, 'measureName',  '',     @ischar);
addParameter(p, 'seed',         42,     @isnumeric);
addParameter(p, 'verbose',      true,   @islogical);
parse(p, dataA, dataB, tvec, varargin{:});
opt = p.Results;
rng(opt.seed);

[nA, nTimeA] = size(dataA);
[nB, nTimeB] = size(dataB);
nTime = numel(tvec);
assert(nTimeA == nTime && nTimeB == nTime, 'time dimension mismatch with tvec.');
nTotal = nA + nB;
pooled = [dataA; dataB];

if opt.verbose
    fprintf('\n--- groupSeparationTest: %s (%s) ---\n', opt.measureName, ternary(opt.paired,'paired','unpaired'));
    fprintf('  nA=%d nB=%d nTime=%d nPerms=%d tail=%s alpha=%.3f corr=%s minG=%.2f\n', ...
        nA, nB, nTime, opt.nPerms, opt.tail, opt.alpha, opt.correction, opt.minHedgesG);
end

%% ---- Restrict to timepoints >= tStart -------------------------------------
testIdx = tvec >= opt.tStart;
tvec = tvec(testIdx); dataA = dataA(:,testIdx); dataB = dataB(:,testIdx);
pooled = pooled(:,testIdx); nTime = sum(testIdx);

%% ---- Observed statistic, effect size, permutation null --------------------
meanA = mean(dataA, 1, 'omitnan');
meanB = mean(dataB, 1, 'omitnan');
meanDiff = meanB - meanA;

permDiffs = zeros(opt.nPerms, nTime);
if opt.paired
    assert(nA == nB, 'groupSeparationTest:paired', 'paired mode needs equal, aligned rows (same sessions).');
    d = dataB - dataA;                              % [n x nTime] paired difference
    sdd = std(d, 0, 1, 'omitnan');
    sdMed = median(sdd, 'omitnan'); sdFloor = 0.01 * sdMed;
    sdd(sdd <= sdFloor | ~isfinite(sdd)) = NaN;
    Jp = 1 - 3 / (4*(nA-1) - 1);                    % small-sample correction, df = n-1
    effectSize = Jp .* (meanDiff ./ sdd);
    for perm = 1:opt.nPerms
        sgn = ones(nA,1); sgn(rand(nA,1) < 0.5) = -1;
        permDiffs(perm,:) = mean(d .* sgn, 1, 'omitnan');
    end
else
    varA = var(dataA, 0, 1, 'omitnan'); varB = var(dataB, 0, 1, 'omitnan');
    pooledSD = sqrt(((nA-1).*varA + (nB-1).*varB) ./ (nTotal - 2));
    sdMed = median(pooledSD, 'omitnan'); sdFloor = 0.01 * sdMed;
    pooledSD(pooledSD <= sdFloor | ~isfinite(pooledSD)) = NaN;
    df = nTotal - 2; J = exp(gammaln(df/2) - log(sqrt(df/2)) - gammaln((df-1)/2));
    effectSize = J .* (meanDiff ./ pooledSD);
    for perm = 1:opt.nPerms
        idx = randperm(nTotal);
        pA = pooled(idx(1:nA), :); pB = pooled(idx(nA+1:end), :);
        permDiffs(perm,:) = mean(pB, 1, 'omitnan') - mean(pA, 1, 'omitnan');
    end
end

%% ---- Raw p-values + BH-FDR ------------------------------------------------
pvals = NaN(1, nTime);
for t = 1:nTime
    obs = meanDiff(t); dist = permDiffs(:, t);
    switch opt.tail
        case 'right'; pvals(t) = mean(dist >= obs);
        case 'left';  pvals(t) = mean(dist <= obs);
        case 'both';  pvals(t) = mean(abs(dist) >= abs(obs));
    end
end
pvals = max(pvals, 1 / opt.nPerms);
qvals = bhFDR(pvals);

%% ---- Cluster correction ---------------------------------------------------
switch opt.tail
    case 'right'; signedObs = meanDiff;      signedPerm = permDiffs;
    case 'left';  signedObs = -meanDiff;     signedPerm = -permDiffs;
    case 'both';  signedObs = abs(meanDiff); signedPerm = abs(permDiffs);
end
switch opt.tail
    case {'right','both'}; primThreshVec = quantile(signedPerm, 1 - opt.primaryAlpha, 1);
    case 'left';           primThreshVec = quantile(signedPerm, opt.primaryAlpha, 1);
end
aboveObs  = signedObs  > primThreshVec;
abovePerm = signedPerm > primThreshVec;

nullMaxMass = zeros(opt.nPerms, 1);
for perm = 1:opt.nPerms
    nullMaxMass(perm) = maxClusterMass(abovePerm(perm,:), signedPerm(perm,:));
end
clusterThresh = quantile(nullMaxMass, 1 - opt.alpha);

clusterMask = false(1, nTime); clusterMass = []; clusterP = [];
if any(aboveObs)
    dd = diff([false aboveObs false]); starts = find(dd==1); ends = find(dd==-1) - 1;
    for c = 1:numel(starts)
        idx = starts(c):ends(c); mass = sum(signedObs(idx));
        pC = max(mean(nullMaxMass >= mass), 1/opt.nPerms);
        clusterMass(end+1) = mass; clusterP(end+1) = pC; %#ok<AGROW>
        if mass > clusterThresh, clusterMask(idx) = true; end
    end
end
if isempty(clusterP), clusterPmin = NaN; else, clusterPmin = min(clusterP); end

%% ---- Active correction -> masks -> onset ----------------------------------
switch opt.correction
    case 'bh';      baseMask = qvals < opt.alpha;
    case 'cluster'; baseMask = clusterMask;
end
sigMask     = baseMask & abs(effectSize) >= opt.minHedgesG;
sigMaskSust = sustainedMask(sigMask, opt.minSustained);
onsetIdx = find(sigMaskSust, 1, 'first');
onsetTime = NaN;
if ~isempty(onsetIdx), onsetTime = tvec(onsetIdx) - opt.onsetOffsetSec; end

%% ---- Output ---------------------------------------------------------------
results = struct('pvals',pvals,'qvals',qvals,'clusterMask',clusterMask, ...
    'clusterMass',clusterMass,'clusterP',clusterP,'clusterPmin',clusterPmin, ...
    'sigMask',sigMask,'sigMaskSust',sigMaskSust,'onsetTime',onsetTime, ...
    'hedgesG',effectSize,'meanDiff',meanDiff,'meanA',meanA,'meanB',meanB, ...
    'tvec',tvec,'measureName',opt.measureName,'params',opt);

if opt.verbose
    fprintf('  Significant (%s + |G|>=%.2f): %d / %d; sustained: %d\n', ...
        opt.correction, opt.minHedgesG, sum(sigMask), nTime, sum(sigMaskSust));
    if ~isnan(onsetTime), fprintf('  --> Onset: %.3f s (cluster p=%.4g)\n', onsetTime, clusterPmin);
    else, fprintf('  --> No sustained separation.\n'); end
end
end

%% =========================================================================
function o = ternary(c,a,b), if c, o=a; else, o=b; end, end
function mass = maxClusterMass(mask, stat)
    dd = diff([false mask(:)' false]); starts = find(dd==1); ends = find(dd==-1) - 1;
    mass = 0; for c = 1:numel(starts), m = sum(stat(starts(c):ends(c))); if m > mass, mass = m; end, end
end
function qvals = bhFDR(pvals)
    m = numel(pvals); [pSorted, sortIdx] = sort(pvals(:)); qSorted = pSorted .* m ./ (1:m)';
    for i = m-1:-1:1, qSorted(i) = min(qSorted(i), qSorted(i+1)); end
    qSorted = min(qSorted, 1); qvals = NaN(size(pvals)); qvals(sortIdx) = qSorted; qvals = reshape(qvals, size(pvals));
end
function mask = sustainedMask(sigIn, minRun)
    mask = false(size(sigIn)); n = numel(sigIn); i = 1;
    while i <= n
        if sigIn(i)
            j = i; while j <= n && sigIn(j), j = j + 1; end
            if j - i >= minRun, mask(i:j-1) = true; end
            i = j;
        else, i = i + 1; end
    end
end
