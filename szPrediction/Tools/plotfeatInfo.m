function [] = plotfeatInfo(X, varargin)
% Plots a matrix with numeric values.

clim = [min(X(:)) max(X(:))];
fmt = '%.3f';

imagesc(X, clim);
colormap(autumn);
ylabel('Patient')
ax.XTickLabel = {'Amp','LL','Alpha','Beta','Delta','Gamma','Theta' 'fxskew' 'fxkurt' 'PSDcorr' 'fycorr' 'fycov' 'HFD'};

% feats = [Amp LL nlEng alpha beta delta gamma theta fxskew fxkurt psdCorr fycorr fycov HFD];

for i = 1:size(X,1)
  for j = 1:size(X,2)
      text(j,i, sprintf(fmt, X(i,j)), ...
           'HorizontalAlignment', 'center', 'color', 'black');    
  end
end
