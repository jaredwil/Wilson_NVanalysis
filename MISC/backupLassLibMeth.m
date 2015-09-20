%Lasso Lib -- just saving this incase I fucked something up

numFolds = 3;
cvIdx    = crossvalind('Kfold', size(retrainLabels,1), numFolds);
numLam   = 100;
lambda   = linspace(10,1e4,numLam);

tLasso_coor = zeros(numFolds,numLam);
tLasso_MSE  = zeros(numFolds,numLam);

for cvIter = 1:numFolds
    disp(['Cross-Validation Progress: ' num2str(cvIter) '/' num2str(numFolds)])
    trainFeatsCV = retrainFeats(cvIdx ~= cvIter,:);
    trainLabelsCV = retrainLabels(cvIdx ~= cvIter,:);
    
    evalFeats = retrainFeats(cvIdx == cvIter,:);
    evalLabels = retrainLabels(cvIdx == cvIter,:);
    %%
    % Train Lasso Regression UseParallel
    w        = cell(1,numLam);
    numIter  = cell(numLam,1);

    tic;
    disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
    parfor_progress(numLam);
    parfor lam = 1:length(lambda)

    %     disp(['Progress: '  num2str(lam) '/' num2str(length(lambda))]);
        [w{lam}, ~, numIter{lam}] = LassoBlockCoordinate(trainFeatsCV,trainLabelsCV,lambda(lam),'maxIter',50000);
    %     [w{lam}, wp, numIter{lam}] = LassoIteratedRidge(trainFeats,trainLabels,lambda(lam),'maxIter',50000);
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
disp(['DONE Training Lasso Model on Patient: ' pt{i}])
lassoTime(i) = toc;
stdMSE = std(tLasso_MSE);

tLasso_MSE = mean(tLasso_MSE,1);

tLasso_coor = mean(tLasso_coor,1);

figure(101)
% plot(lambda, tLasso_MSE,'r.-')
errorbar(lambda, tLasso_MSE,stdMSE,'r.-')

% minIdx = tINFO.IndexMinMSE;
% seIdx  = tINFO.Index1SE;
corrIdx = tLasso_coor == max(tLasso_coor);
minIdx = tLasso_MSE == min(tLasso_MSE);
testIdx = (tLasso_MSE - stdMSE) == min(tLasso_MSE - stdMSE);

figure(102)
title('Number of Non-Zero Features vs Resulting Correlation')
h = plot(lambda,tLasso_coor*100,'r.-');
% xlabel('Number of Feature Weights Zeroed in Lasso Model')
xlabel('Lambda')
ylabel('Resulting Test Correlation (%)')
grid on;
hold on;
h = plot(lambda(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
h = vline(lambda(minIdx),'g:');
h = plot(lambda(corrIdx),tLasso_coor(corrIdx)*100,'bs','LineWidth',4);
h = plot(lambda(testIdx),tLasso_coor(testIdx)*100,'ms','LineWidth',4);
% h = vline(dzT(seIdx),'b:');
plotName = ['H:\jaredwil\Lasso Results\' pt{i} '_retrainCorrLib'];
saveas(h,plotName,'jpg')

[bestLasso_test, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(testIdx),'maxIter',50000);
[bestLasso_corr, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(corrIdx),'maxIter',50000);
[bestLasso_Min, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(minIdx),'maxIter',50000);


numFeats_test = sum(bestLasso_test ~= 0);
numFeats_corr = sum(bestLasso_corr ~= 0);
numFeats_Min = sum(bestLasso_Min ~= 0);

bestInt_corr = mean(retrainLabels) -  mean(retrainFeats*bestLasso_test);
bestInt_test = mean(retrainLabels) -  mean(retrainFeats*bestLasso_corr);
bestInt_Min  = mean(retrainLabels) -  mean(retrainFeats*bestLasso_Min);

