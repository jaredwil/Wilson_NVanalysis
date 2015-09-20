function [ f, int, numFeats, solTime ] = trainLassoReg( feats,labels,numFolds,numlambda, modCh )
%Train a Lasso Model using external Lasso Library using CV
%   Detailed explanation goes here

modCh = lower(modCh); 
%this choses which model to return
%   min - min MSE
%   corr - max correlation
%   se   - 1 standard deviation away from min MSE

%Reserve some labels for evaluation set
tr = length(labels)*0.7;
evalFeats = feats(tr+1:end,:);
evalLabels = labels(tr+1:end);
feats = feats(1:tr,:);
labels = labels(1:tr); 

tic;
opts = statset('UseParallel',true);  %Do dis shit in parallel!!!  use number of availble workers for CV
[fLasso, tINFO] = lasso(feats, labels,'CV',numFolds,'NumLambda',numlambda,'Options',opts);
solTime = toc;

dzT = tINFO.DF; 
tInt = tINFO.Intercept;
% xAlpha = xINFO.Alpha;
tInt2 = repmat(tInt, size(feats,1),1);
utLasso = feats*fLasso + tInt2;
tTrain = repmat(labels,1,100);
tLasso_coor = corr(tTrain,utLasso);
tLasso_coor = tLasso_coor(1,:);

figure(11)
plot(dzT,tLasso_coor*100,'r.');
xlabel('Number of Non-Zero Feature Weights in Lasso Model')
ylabel('Resulting Test Correlation (%)')
title('TRAIN COULD BE MISLEADING')
grid on;
hold on;

tInt2 = repmat(tInt, size(evalFeats,1),1);
utLasso = evalFeats*fLasso + tInt2;
tTrain = repmat(evalLabels,1,100);
tLasso_coor = corr(tTrain,utLasso);
tLasso_coor = tLasso_coor(1,:);

minIdx = tINFO.IndexMinMSE;
seIdx  = tINFO.Index1SE;
corrIdx = tLasso_coor == max(tLasso_coor);

figure(12)
plot(dzT,tLasso_coor*100,'r.');
xlabel('Number of Non-Zero Feature Weights in Lasso Model')
ylabel('Resulting Test Correlation (%)')
title('EVALUATION')
grid on;
hold on;
h = plot(dzT(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
h = vline(dzT(minIdx),'g:');
h = plot(dzT(seIdx),tLasso_coor(seIdx)*100,'bs','LineWidth',4);
h = vline(dzT(seIdx),'b:');

figure(14)
h2 = lassoPlot(fLasso,tINFO,'PlotType','CV');
% Use a log scale for MSE to see small MSE values better
% set(gca,'YScale','log');
% plotName = ['H:\jaredwil\Lasso Results\matLasso_history\' pt{i} '_mseRes_CCS'];

bestLasso_corr = fLasso(:,corrIdx);
bestLasso_1SE = fLasso(:,seIdx);
bestLasso_Min = fLasso(:,minIdx);

int_corr = tInt(corrIdx); 
int_1SE = tInt(seIdx); 
int_Min = tInt(minIdx); 

numFeats_corr = dzT(corrIdx); 
numFeats_1SE = dzT(seIdx); 
numFeats_Min = dzT(minIdx);

switch modCh
    case 'corr'
        f = bestLasso_corr;
        int = int_corr;
        numFeats = numFeats_corr;
    case 'min'
        f = bestLasso_Min;
        int = int_Min;
        numFeats = numFeats_Min;        
    case 'se'
        f = bestLasso_1SE;
        int = int_1SE;
        numFeats = numFeats_1SE;       
end


end

