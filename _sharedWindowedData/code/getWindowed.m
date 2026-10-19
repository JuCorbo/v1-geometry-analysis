function [X, info] = getWindowed(S, scheme)
%GETWINDOWED  The single accessor for the shared windowed intermediate.
%
%   Every script obtains its data through here, so the windowing, the baseline
%   subtraction and the sigma floor are defined in exactly one place.
%
%   S       a loaded session struct (from load(<...>_win.mat)) or a path to one
%   scheme  'raw'       raw window sums, no normalisation           [Fig 4 capacity]
%           'sub'       per-trial baseline subtracted, no z-score
%           'baseline'  baseline z-score, sigma floored at 0.12     [Fig 1, Fig 3-native, Fig 2 supp, Fig 4]
%           'global'    baseline subtracted then GLOBAL z-score     [Fig 2 main]
%           'baseline_nofloor'  as 'baseline' but with the unfloored sigma
%                               (for the with/without-floor sensitivity analysis)
%
%   Sigma floor = 0.12 (see the Methods and the Figure1/Supp_sigma_floor
%   diagnostic). The floor is applied here, to sd_base_unfloored; the cache's
%   stored sd_base field is not used.
%
%   X     double [nTrials x nCells x nWin]
%   info  struct with tvec, tvec_centre, bins, trial/cell metadata, provenance
%
%   Notes
%     * every bin is a full 3 frames: (c-2):c for c = 9:3:84  -> {7,8,9} .. {82,83,84}
%     * because all bins are the same width, baseline-subtracted = X - 3*b
%     * tvec is the LAST-frame time; tvec_centre is the true centre of a
%       backward window (centers+(winFrames-1)/2 would be 2 frames late).

    if ischar(S) || isstring(S), S = load(S); end
    if nargin < 2, scheme = 'baseline'; end

    Xr = double(S.X);
    nW = size(Xr,3);

    switch lower(scheme)
        case 'raw'
            X = Xr;
        case 'sub'
            X = Xr - 3*repmat(double(S.b),1,1,nW);
        case 'baseline'
            % baseline z-score, sigma floored at 0.12, from the unfloored SD.
            sd = double(S.sd_base_unfloored); sd(~isfinite(sd) | sd < 0.12) = 0.12;
            X = zscoreBy(Xr, S.mu_base, sd);
        case 'baseline_nofloor'
            sd = S.sd_base_unfloored; sd(~isfinite(sd) | sd<=0) = 1;   % avoid /0 only
            X = zscoreBy(Xr, S.mu_base, sd);
        case 'global'
            X = zscoreBy(Xr - 3*repmat(double(S.b),1,1,nW), S.mu_glob_sub, S.sd_glob_sub);
        otherwise
            error('getWindowed:scheme','unknown scheme "%s"', scheme);
    end

    if nargout > 1
        info = struct('tvec',S.tvec,'tvec_centre',S.tvec_centre,'bins',{S.bins}, ...
            'centers',S.centers,'baselineBins',S.baselineBins, ...
            'trialID',S.trialID,'Visual_Stim',S.Visual_Stim,'Response',S.Response, ...
            'Outcome',{S.Outcome},'cellID',S.cellID,'OSI',S.OSI,'PrefOri',S.PrefOri, ...
            'session',S.session,'day',S.day,'group',S.group,'DB',S.DB, ...
            'dt',S.dt,'winFrames',S.winFrames,'sd_floor',S.sd_floor, ...
            'provenance',S.provenance);
    end
end

function Z = zscoreBy(X, mu, sd)
    nW = size(X,3);
    Z = (X - repmat(reshape(mu,1,[],1),size(X,1),1,nW)) ./ ...
             repmat(reshape(sd,1,[],1),size(X,1),1,nW);
end
