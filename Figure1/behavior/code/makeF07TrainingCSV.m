%% makeF07TrainingCSV.m — stage session-resolved training behavior (Zenodo) into the deliverable
%  Source: the 2025-paper code deposit DataFiles (session-resolved d'/hit/FA per animal,
%          Beh rig = pre-imaging training), C:\...\MATLAB\Zenodo\DataFiles.
%  Writes (SQL-free, self-contained plot data):
%    data/fig1_c_training_exampleF07.csv        — the Figure-1c example (F07), Beh phase
%    data/trainingSession_allAnimals_Beh.csv     — all 10 trained mice, Beh phase (for easy swaps)
%  Session-resolved (one d' per training session), unlike the trial-resolved rolling curve that
%  only B10/C02/C05 have (raw NAS tally files); the F-cohort exists only at session resolution.
clear; clc;
codeDir=fileparts(mfilename('fullpath')); if isempty(codeDir), codeDir=pwd; end
aDir=fileparts(codeDir); dataDir=fullfile(aDir,'data'); if ~exist(dataDir,'dir'), mkdir(dataDir); end
% BUILDER INPUT (not distributed): local copy of the public 2025-paper code
% deposit's DataFiles folder (Zenodo, record 8109858 code deposit). The CSVs in
% ../data are this script's product; set zdir to your download to rebuild.
if ~exist('zdir','var'), error('Set zdir to your local copy of the deposit DataFiles folder (see header).'); end

ld=@(f)load(fullfile(zdir,f)); g=@(S)S.(char(fieldnames(S)));
DB=g(ld('dprimeBeh.mat')); HB=g(ld('hitRateBeh.mat')); FB=g(ld('faRateBeh.mat'));
pick=@(C,nm)C{find(string(C(:,1))==nm,1),2};
toPct=@(v) (100*v(:)').*(max(v(:))<=1) + (v(:)').*(max(v(:))>1);   % fractions->%, pass-through if already %

names=string(DB(:,1));
EX="F07";

% --- example CSV ---
d=pick(DB,EX); d=d(:)'; h=toPct(pick(HB,EX)); f=toPct(pick(FB,EX)); s=(1:numel(d))';
Tex=table(s, h(:), f(:), d(:), 'VariableNames',{'session','hitPct','faPct','dprime'});
writetable(Tex, fullfile(dataDir,'fig1_c_training_exampleF07.csv'));
fprintf('wrote fig1_c_training_exampleF07.csv (%s, %d sessions, d'' %.2f->%.2f)\n', EX, numel(d), d(1), d(end));

% --- all-animals CSV (Beh phase) ---
rows=table();
for i=1:numel(names)
  nm=names(i); dd=pick(DB,nm); dd=dd(:)'; hh=toPct(pick(HB,nm)); ff=toPct(pick(FB,nm));
  n=numel(dd);
  rows=[rows; table(repmat(nm,n,1),(1:n)',hh(:),ff(:),dd(:), ...
    'VariableNames',{'animal','session','hitPct','faPct','dprime'})]; %#ok<AGROW>
end
writetable(rows, fullfile(dataDir,'trainingSession_allAnimals_Beh.csv'));
fprintf('wrote trainingSession_allAnimals_Beh.csv (%d animals, %d session-rows)\n', numel(names), height(rows));
