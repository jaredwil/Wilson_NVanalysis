function [ f, int, numFeats, solTime ] = trainLassoLib( feats,labels,numFolds,lambda, modCh )
%Train a Lasso Model using external Lasso Library using CV
%   Detailed explanation goes here

modCh = lower(modCh);

cvIdx    = crossvalind('Kfold', size(labels,1), numFolds);
numLam   = length(lambda);

tLasso_coor = zeros(numFolds,numLam);
tLasso_MSE  = zeros(numFolds,numLam);

for cvIter = 1:numFolds
    disp(['Cross-Validation Progress: ' num2str(cvIter) '/' num2str(numFolds)])
    trainFeatsCV = feats(cvIdx ~= cvIter,:);
    trainLabelsCV = labels(cvIdx ~= cvIter,:);
    
    evalFeats = feats(cvIdx == cvIter,:);
    evalLabels = labels(cvIdx == cvIter,:);
    %%
    % Train Lasso Regression UseParallel
    w        = cell(1,numLam);
    numIter  = cell(numLam,1);

    tic;
    parfor_progress(numLam);
    parfor lam = 1:length(lambda)
        %change this functiion manually if you wish to use a dif method
        [w{lam}, ~, numIter{lam}] = LassoBlockCoordinate(trainFeatsCV,trainLabelsCV,lambda(lam),'maxIter',50000);
        parfor_progress;
    end
    parfor_progress(0);

    w = cell2mat(w);

    dzT = sum(w == 0,1);
    %compute intercepts
    tInt = repmat(mean(trainLabelsCV),1,numLam) - mean(trainFeatsCV*w,1);

    tInt2 = repmat(tInt, size(evalFeats,1),1);

    utLasso = evalFeats*w + tInt2;

    tTrain = repmat(evalLabels,1,numLam);

    tmp = corr(tTrain,utLasso);
    tLasso_coor(cvIter,:) = tmp(1,:);

    tLasso_MSE(cvIter,:) = mean(bsxfun(@minus,utLasso,evalLabels).^2,1);

end
%BOOM!!! Done.
solTime = toc;

stdMSE = std(tLasso_MSE);
tLasso_MSE = mean(tLasso_MSE,1);

tLasso_coor = mean(tLasso_coor,1);

figure(101)
errorbar(lambda, tLasso_MSE,stdMSE,'r.-')

% minIdx = tINFO.IndexMinMSE;
% seIdx  = tINFO.Index1SE;
corrIdx = tLasso_coor == max(tLasso_coor);
minIdx = tLasso_MSE == min(tLasso_MSE);

figure(102)
title('Number of Non-Zero Features vs Resulting Correlation')
h = plot(lambda,tLasso_coor*100,'r.-');
% xlabel('Number of Feature Weights Zeroed in Lasso Model')
xlabel('Lambda')
ylabel('Resulting Test Correlation (%)')
grid on;
hold on;
h = plot(dzT(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
h = vline(dzT(minIdx),'g:');
% h = plot(dzT(seIdx),tLasso_coor(seIdx)*100,'bs','LineWidth',4);
% h = vline(dzT(seIdx),'b:');
% plotName = ['H:\jaredwil\Lasso Results\' pt{i} '_corrRes'];
% saveas(h,plotName,'jpg')

testIdx = 53;

[bestLasso_test, ~, ~] = LassoBlockCoordinate(feats,labels,lambda(testIdx),'maxIter',50000);
[bestLasso_corr, ~, ~] = LassoBlockCoordinate(feats,labels,lambda(corrIdx),'maxIter',50000);
[bestLasso_Min, ~, ~] = LassoBlockCoordinate(feats,labels,lambda(minIdx),'maxIter',50000);


numFeats_test = sum(bestLasso_test ~= 0);
numFeats_corr = sum(bestLasso_corr ~= 0);
numFeats_Min = sum(bestLasso_Min ~= 0);

int_test = mean(labels) -  mean(feats*bestLasso_test);
int_corr = mean(labels) -  mean(feats*bestLasso_corr);
int_Min  = mean(labels) -  mean(feats*bestLasso_Min);

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
        f = bestLasso_test;
        int = int_test;
        numFeats = numFeats_test;       
end


end

