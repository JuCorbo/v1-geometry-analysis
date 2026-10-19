%% kNN short-circuit quantification: global-z vs baseline-z neighbour graphs.
% Supplementary-methods point: under a GLOBAL z-score the
% per-neuron variance weighting self-attenuates discriminative neurons, so nearest-
% neighbour graphs acquire cross-cluster edges (Isomap/UMATO short-circuit); a
% BASELINE z-score preserves cluster separability. Both normalizations come from
% the shared cache via getWindowed('global'|'baseline').
clear; clc;
codeDir=fileparts(mfilename('fullpath'));
figDir=fileparts(fileparts(fileparts(codeDir)));           % Figures_final
sharedCode=fullfile(figDir,'_sharedWindowedData','code');
sharedData=fullfile(figDir,'_sharedWindowedData','data');
addpath(sharedCode);

grp={'naive','Naive';'trained','Trained'}; days=[1 3 5]; k=12; goOri=45;
fprintf('%-28s | GLOBAL-z         | BASELINE-z\n','session (day)');
fprintf('%-28s | cross%%   b/w    | cross%%   b/w\n','');
for gi=1:2
  for d=days
    ff=dir(fullfile(sharedData,grp{gi,1},sprintf('D%d',d),'*_win.mat'));
    if isempty(ff), continue; end
    S=load(fullfile(ff(1).folder,ff(1).name));
    [Xg,info]=getWindowed(S,'global'); [Xb,~]=getWindowed(S,'baseline');   % [nTr x nC x nW]
    vs=double(info.Visual_Stim(:)); tv=double(info.tvec(:));
    sv=unique(vs(isfinite(vs))); go=goOri; if ~any(sv==goOri), go=sv(1); end; no=sv(sv~=go); no=no(1);
    m=(vs==go|vs==no);
    evW=find(tv>0.2 & tv<3.0); [~,swi]=min(abs(tv-1.5));
    % event-window point cloud (trials x ev-windows) and 1.5 s snapshot
    Xg_ev=[]; Xb_ev=[]; L_ev=[];
    for w=evW', Xg_ev=[Xg_ev; Xg(m,:,w)]; Xb_ev=[Xb_ev; Xb(m,:,w)]; L_ev=[L_ev; vs(m)]; end %#ok<AGROW>
    Xg_sm=Xg(m,:,swi); Xb_sm=Xb(m,:,swi); L_sm=vs(m);
    out=zeros(2,2);
    for which=1:2
      if which==1, X=Xg_ev; Xs=Xg_sm; else, X=Xb_ev; Xs=Xb_sm; end
      [idx,~]=knnsearch(X,X,'K',k+1); idx=idx(:,2:end);
      crossf=100*mean(mean(L_ev(idx)~=L_ev,2));
      mg=mean(Xs(L_sm==go,:),1); mn=mean(Xs(L_sm==no,:),1);
      wthn=mean([sqrt(mean(sum((Xs(L_sm==go,:)-mg).^2,2))),sqrt(mean(sum((Xs(L_sm==no,:)-mn).^2,2)))]);
      out(which,:)=[crossf, norm(mg-mn)/wthn];
    end
    nm=char(info.session);
    fprintf('%-28s | %5.1f%%  %5.2f  | %5.1f%%  %5.2f\n', sprintf('%s (D%d)',nm(1:min(20,end)),d), out(1,1),out(1,2),out(2,1),out(2,2));
  end
end
