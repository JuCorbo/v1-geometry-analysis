function T = behavDepositTrials(group, visualOnly)
%BEHAVDEPOSITTRIALS  Per-trial behavioral scalars from the public deposit (SQL-free).
%   T = behavDepositTrials()                -> Trained group, Visual block only
%   T = behavDepositTrials(group)           -> 'Trained' (default) or 'Naive'
%   T = behavDepositTrials(group, false)    -> keep all blocks (incl. Orientation Tuning)
%
%   Returns a table with columns:
%     Trial (double), Experiment (string), Behav_Cond (string),
%     Block (string), Visual_Stim (double), Response (double 0/1), Outcome (string)
%
%   Source: _sharedBehavior/data/<group>_behav_slim.csv, a scalar-column
%   slice of <group>_TrialTable.csv from DatasetsCorbo2024.rar (Zenodo deposit).
%   The deposit's truncated Lick/Loco/Pupil trace blobs are NOT included here;
%   only the per-trial scalars needed for lick-yes/no and d'/FA, which are exact.

    if nargin < 1 || isempty(group),      group = 'Trained'; end
    if nargin < 2 || isempty(visualOnly), visualOnly = true;  end
    group = char(group);

    here = fileparts(mfilename('fullpath'));
    csv  = fullfile(here, '..', 'data', [group '_behav_slim.csv']);
    if ~isfile(csv)
        error('behavDepositTrials:missing', ...
            'Slim behavioral table not found:\n  %s\nExtract it from DatasetsCorbo2024.rar.', csv);
    end

    T = readtable(csv, 'TextType','string', 'Delimiter',',');
    T.Experiment = string(T.Experiment);
    T.Behav_Cond = string(T.Behav_Cond);
    T.Block      = string(T.Block);
    T.Outcome    = string(T.Outcome);

    if visualOnly
        T = T(T.Block == "Visual", :);
    end
end
