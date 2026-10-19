function [dprime, pFAall, pHitall, sessAllDays, days] = behavDprimeBySession(group, goStim, nSessCols)
%BEHAVDPRIMEBYSESSION  Per-session d' / FA from the deposit CSV (SQL-free).
%
%   [dprime, pFAall, pHitall, sessAllDays, days] = behavDprimeBySession(group, goStim, nSessCols)
%     group     : 'Trained' (default) or 'Naive'
%     goStim    : Go orientation (default 45)
%     nSessCols : width of the output matrices (default 10; days with fewer
%                 sessions are NaN-padded)
%
%   Outputs:
%     dprime, pFAall, pHitall : [6 x nSessCols], row = day D1..D6, col = session index
%     sessAllDays             : {6 x 1} cell, each a cellstr of Experiment IDs (sorted)
%     days                    : ["D1".."D6"]
%
%   d' = norminv(clip(pHit)) - norminv(clip(pFA)), clip to [0.01, 0.99];
%   pHit over Go trials, pFA over NoGo (all non-Go visual stimuli). Visual block only.

    if nargin < 1 || isempty(group),     group     = 'Trained'; end
    if nargin < 2 || isempty(goStim),    goStim    = 45;        end
    if nargin < 3 || isempty(nSessCols), nSessCols = 10;        end

    T    = behavDepositTrials(group, true);      % Visual block only
    days = "D" + string(1:6);

    dprime      = nan(6, nSessCols);
    pFAall      = nan(6, nSessCols);
    pHitall     = nan(6, nSessCols);
    sessAllDays = cell(6, 1);

    clip = @(p) max(min(p, 0.99), 0.01);

    for di = 1:6
        Td = T(T.Behav_Cond == days(di), :);
        if isempty(Td), sessAllDays{di} = {}; continue; end

        sess = unique(Td.Experiment);            % sorted
        sessAllDays{di} = cellstr(sess);
        if numel(sess) > nSessCols
            warning('behavDprimeBySession:width', ...
                'D%d has %d sessions > nSessCols=%d; widen nSessCols.', di, numel(sess), nSessCols);
        end

        for si = 1:min(numel(sess), nSessCols)
            Te   = Td(Td.Experiment == sess(si), :);
            go   = Te.Visual_Stim == goStim;
            nogo = ~go;                          % all non-Go visual stimuli
            pHit = mean(Te.Response(go),   'omitnan');
            pFA  = mean(Te.Response(nogo), 'omitnan');

            pHitall(di, si) = pHit;
            pFAall(di, si)  = pFA;
            dprime(di, si)  = norminv(clip(pHit)) - norminv(clip(pFA));
        end
    end
end
