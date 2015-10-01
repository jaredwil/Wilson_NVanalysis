function [] = plotchannelInfo(X, varargin)
% Plots a matrix with numeric values.

clim = [min(X(:)) max(X(:))];
fmt = '%.2f';

r = size(X,1);
c = size(X,2);

imagesc(X, clim);
colormap(autumn);
ylabel('Patient')
xlabel('Channel')
set(gca,'XTick',linspace(1,c,c),'YTick',linspace(1,r,r));
ax = gca;
% ax.XTickLabel = {'Amp','LL','NLenergy','Alpha','Beta','Delta','Gamma','Theta' 'fxskew' 'fxkurt' 'PSDcorr' 'fycorr' 'fycov' 'HFD'};

% feats = [Amp LL nlEng alpha beta delta gamma theta fxskew fxkurt psdCorr fycorr fycov HFD];

for i = 1:size(X,1)
  for j = 1:size(X,2)
      text(j,i, sprintf(fmt, X(i,j)), ...
           'HorizontalAlignment', 'center', 'color', 'black');    
  end
end
