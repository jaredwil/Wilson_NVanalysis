%Jared Wilson 
%Retrain Sz based on Correlation
%Find the correlation of the predicted time to sz. for each seizure and see
%if retraining increases the correlation.

%%
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\30secFeats'))  %this is where .mat file are contained on local comp
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\5minFeats'))  %this is where .mat file are contained on local comp
addpath(genpath('H:\jaredwil\szPred_feats\Win_5min'))  %this is where .mat file are contained on local comp


%%
% Define algorithm specifics 
usernm = 'jaredwil';
pswdBin = 'jar_ieeglogin.bin';
trPct = 0.7;
winLen = 30;
winDisp = 30;
szHorizon = 2; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%%%%%%%%%%%%%%%%%%%%%% start par pool%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Establish parpool
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)        %if not p will contain all info about current pool;
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
%initialize paralell pool if available
if(poolsize == 0)

    myCluster = parcluster('local');
    numWork = myCluster.NumWorkers;

    if(numWork <= 2)  %processing is done on a laptop so don't do it in parallel
        parpool('local',1)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;
    elseif(numWork > 10) %limit the number of workers to 6
        parpool('local',10)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;
    else  %set up a parallel pool with max number of workers available between 2 and 6
        parpool(myCluster)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lassoTime = zeros(1,numel(pt));
%%
%begin function
%loop through all pts
for i = 7%:numel(pt);  %%%%%TEMPORARY for debug
       
%     label = [pt{i} '_szPred_30secFeats.mat'];
    label = [pt{i} '_szPred_5minFeats.mat'];

    tylabel = [pt{i} '_szTypeLabels5.mat'];

    try   
        load(label);
        load(tylabel);
    catch
        disp('This pt. does not have a save .mat to load from')
        continue;

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% get feature data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    train   = data.train;
    test    = data.test;
    % avgFeats = data.mean;
    % stdFeats = data.std;

    %isolate labels and feats in training data
    trainLabels = train(:,1);
    trainFeats = train(:,2:end);
    testLabels = test(:,1);
    testFeats = test(:,2:end);

    %Remove intericatal Data from training data & Testing Data
    trainFeats(trainLabels > szHorizon*60*60,:) = [];
    trainLabels(trainLabels > szHorizon*60*60,:) = [];
    testFeats(testLabels > szHorizon*60*60,:) = [];
    testLabels(testLabels > szHorizon*60*60,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%% get sz time and type labels%%%%%%%%%%%%%%%%%%%%

    trainSzT = szType.labels.train;
    testSzT = szType.labels.test;
    
    szTimes = [szType.szT.train;szType.szT.test];
    
    szTimes(:,1) = szTimes(:,1)/60/60/24;
    
    avgFeats = mean(trainFeats,1);
    stdFeats = std(trainFeats,[],1);

    allFeat = [trainFeats; testFeats];
    allFeat = bsxfun(@rdivide, bsxfun(@minus,allFeat,avgFeats), stdFeats);

    allLabs  = [trainLabels;testLabels];
    allSzLab  = [trainSzT;testSzT];
    
    
%%%%%%%%%%%%%% REMOVE SOME SZ %%%%%%%%%%%%%%%%%%%%%%%%%%%
%ONLY LOOK AT CCS
allFeat(allSzLab(:,2)~=1,:) = [];
allLabs(allSzLab(:,2)~=1) = [];
allSzLab(allSzLab(:,2)~=1,:) = [];

% %%%%%%%%%%%%%%%%%%%%%%%RESHAPE FEATURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1:2; %number of history samples
% samples = size(allFeat,1);
% features = size(allFeat,2);
% startSAMP = max(N)+1;
% M = (samples-length(N));  %number of time bins
% %create 'R' matrix for linear regression algorithm
% r = zeros(M, features*length(N)+1);
% for resIdx = 1:M
%     temp = allFeat(startSAMP + (resIdx-1) - N,:);   %temp is a temporary matrix    
%     r(resIdx,:) = [1 temp(:)'];
% end
% 
% allLabs = allLabs(startSAMP:end,:);
% allSzLab = allSzLab(startSAMP:end,:);
% 
% [~,IdxRemove] = findpeaks(allLabs);
% 
% r(IdxRemove,:) = [];
% allSzLab(IdxRemove,:) = [];
% allLabs(IdxRemove,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%load feature weights
% load('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\resultsV1\NVC1001_23_005_bestLasso_CCS.mat')

% featLab = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Res LassoLib\' pt{i} '_bestLasso5.mat'];
featLab = ['H:\jaredwil\Lasso Results\Res 8-29\' pt{i} '_bestLasso5.mat'];

load(featLab);


f = lassoRes.coef.totmin;
tInt = lassoRes.int(3);
tInt2 = repmat(tInt, size(allFeat,1),1);


predLab = allFeat*f + tInt2;

%visual assesment
figure(1)
plot(allLabs);
hold on;
plot(smooth(predLab,5));


    %average feature across sz do pca then cluster. 
    [numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 
    
    szCorr = zeros(numel(numSz),1);
    timeDay = zeros(numel(numSz),1);
    timeSec = zeros(numel(numSz),1);
    for sz = 1:numel(numSz);
%         szAvg_feat(sz,:) = [mean(allFeat(uIdx == sz,:),1) szTimes(sz,1)];
        szCorr(sz) = corr(allLabs(uIdx == sz),predLab(uIdx == sz));
        timeDay(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(allSzLab(uIdx == sz,1));

    end
        
    %visualize what I am doing
    figure(2)
    plot(timeDay,szCorr*100,'r.','MarkerSize',20)
    hold on;
    grid on;
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1);
    set(gcf,'position',get(0,'screensize'));
    vline(timeDay)
    vline(szType.szT.test(1)/60/60/24,'b')
    hline(60,'m');
    ylabel('Correlation (%)')
    xlabel('Time (days)')
    title('Correlation of Predicted Time to Seizure for Each Seizure')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Original corr by sz has been found now threshold and keep sz's with high
% correlation then retrain

corrTh = 0.6; %60 percent

corrLab = zeros(size(allLabs));
isTrain = zeros(size(allLabs));

for sz = 1:numel(numSz);
    corrLab(uIdx == sz) = szCorr(sz);                                   %repeat correlation for a label
    isTrain(uIdx == sz) = sum((timeSec(sz) == szType.szT.train(:,1)));  %sz was orignally trained on
end
isTrain = logical(isTrain);

idxTrain = (corrLab > corrTh) & isTrain; %these are the feature idx to train a new model on

retrainFeats  = allFeat(idxTrain,:);
retrainLabels = allLabs(idxTrain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method Using Lasso Lib
% numFolds = 5;
% cvIdx    = crossvalind('Kfold', size(retrainLabels,1), numFolds);
% numLam   = 100;
% lambda   = linspace(10,1e6,numLam);
% 
% tLasso_coor = zeros(numFolds,numLam);
% tLasso_MSE  = zeros(numFolds,numLam);
% 
% for cvIter = 1:numFolds
%     disp(['Cross-Validation Progress: ' num2str(cvIter) '/' num2str(numFolds)])
%     trainFeatsCV = retrainFeats(cvIdx ~= cvIter,:);
%     trainLabelsCV = retrainLabels(cvIdx ~= cvIter,:);
%     
%     evalFeats = retrainFeats(cvIdx == cvIter,:);
%     evalLabels = retrainLabels(cvIdx == cvIter,:);
%     %%
%     % Train Lasso Regression UseParallel
%     w        = cell(1,numLam);
%     numIter  = cell(numLam,1);
% 
%     tic;
%     disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
%     parfor_progress(numLam);
%     parfor lam = 1:length(lambda)
% 
%     %     disp(['Progress: '  num2str(lam) '/' num2str(length(lambda))]);
%         [w{lam}, ~, numIter{lam}] = LassoBlockCoordinate(trainFeatsCV,trainLabelsCV,lambda(lam),'maxIter',50000);
%     %     [w{lam}, wp, numIter{lam}] = LassoIteratedRidge(trainFeats,trainLabels,lambda(lam),'maxIter',50000);
%         parfor_progress;
% 
%     end
%     parfor_progress(0);
% 
% 
% 
%     w = cell2mat(w);
% 
%     dzT = sum(w == 0,1);
%     %compute intercepts
%     tInt = repmat(mean(trainLabelsCV),1,numLam) - mean(trainFeatsCV*w,1);
% 
%     tInt2 = repmat(tInt, size(evalFeats,1),1);
% 
%     utLasso = evalFeats*w + tInt2;
% 
%     tTrain = repmat(evalLabels,1,numLam);
% 
%     tmp = corr(tTrain,utLasso);
%     tLasso_coor(cvIter,:) = tmp(1,:);
% 
%     tLasso_MSE(cvIter,:) = mean(bsxfun(@minus,utLasso,evalLabels).^2,1);
% 
% end
% %BOOM!!! Done.
% disp(['DONE Training Lasso Model on Patient: ' pt{i}])
% lassoTime(i) = toc;
% stdMSE = std(tLasso_MSE);
% 
% tLasso_MSE = mean(tLasso_MSE,1);
% 
% tLasso_coor = mean(tLasso_coor,1);
% 
% figure(101)
% plot(lambda, tLasso_MSE,'r.-')
% 
% % minIdx = tINFO.IndexMinMSE;
% % seIdx  = tINFO.Index1SE;
% corrIdx = tLasso_coor == max(tLasso_coor);
% minIdx = tLasso_MSE == min(tLasso_MSE);
% 
% figure(102)
% title('Number of Non-Zero Features vs Resulting Correlation')
% h = plot(lambda,tLasso_coor*100,'r.-');
% % xlabel('Number of Feature Weights Zeroed in Lasso Model')
% xlabel('Lambda')
% ylabel('Resulting Test Correlation (%)')
% grid on;
% hold on;
% h = plot(dzT(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
% h = vline(dzT(minIdx),'g:');
% % h = plot(dzT(seIdx),tLasso_coor(seIdx)*100,'bs','LineWidth',4);
% % h = vline(dzT(seIdx),'b:');
% % plotName = ['H:\jaredwil\Lasso Results\' pt{i} '_corrRes'];
% % saveas(h,plotName,'jpg')
% 
% testIdx = 53;
% 
% [bestLasso_test, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(testIdx),'maxIter',50000);
% [bestLasso_corr, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(corrIdx),'maxIter',50000);
% [bestLasso_Min, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(minIdx),'maxIter',50000);
% 
% 
% % bestLasso_corr = w(:,corrIdx);
% % % bestLasso_1SE = fLasso(:,seIdx);
% % bestLasso_Min = w(:,minIdx);
% 
% numFeats_test = sum(bestLasso_test ~= 0);
% numFeats_corr = sum(bestLasso_corr ~= 0);
% % numFeats_1SE = dzT(seIdx); 
% numFeats_Min = sum(bestLasso_Min ~= 0);
% 
% int_test = mean(retrainLabels) -  mean(retrainFeats*bestLasso_test);
% int_corr = mean(retrainLabels) -  mean(retrainFeats*bestLasso_corr);
% int_Min  = mean(retrainLabels) -  mean(retrainFeats*bestLasso_Min);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %Other method
%split up for a simple 2 fold CV
tr = length(retrainLabels)*0.7;
evalFeats = retrainFeats(tr+1:end,:);
evalLabels = retrainLabels(tr+1:end);
retrainFeats = retrainFeats(1:tr,:);
retrainLabels = retrainLabels(1:tr); 

tic;
disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
opts = statset('UseParallel',true);  %Do dis shit in parallel!!!  use number of availble workers for CV
[fLasso, tINFO] = lasso(retrainFeats, retrainLabels,'CV',5,'Options',opts);
%BOOM!!!
disp(['DONE Training Lasso Model on Patient: ' pt{i}])
dzT = tINFO.DF; 
tInt = tINFO.Intercept;
% xAlpha = xINFO.Alpha;
tInt2 = repmat(tInt, size(retrainFeats,1),1);
utLasso = retrainFeats*fLasso + tInt2;
tTrain = repmat(retrainLabels,1,100);
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

bestInt_corr = tInt(corrIdx); 
bestInt_1SE = tInt(seIdx); 
bestInt_Min = tInt(minIdx); 

numFeats_corr = dzT(corrIdx); 
numFeats_1SE = dzT(seIdx); 
numFeats_Min = dzT(minIdx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%NOW REDO THE FIREST PART WITH THE NEW STUFFS

% 
% f = bestLasso_test;
% tInt = int_test;
f = bestLasso_corr;
tInt = bestInt_corr;
tInt2 = repmat(tInt, size(allFeat,1),1);


predLab = allFeat*f + tInt2;

%visual assesment
figure(100)
plot(allLabs);
hold on;
plot(smooth(predLab,5));


    %average feature across sz do pca then cluster. 
    [numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 
    
    szCorr = zeros(numel(numSz),1);
    timeDay = zeros(numel(numSz),1);
    timeSec = zeros(numel(numSz),1);
    for sz = 1:numel(numSz);
%         szAvg_feat(sz,:) = [mean(allFeat(uIdx == sz,:),1) szTimes(sz,1)];
        szCorr(sz) = corr(allLabs(uIdx == sz),predLab(uIdx == sz));
        timeDay(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(allSzLab(uIdx == sz,1));

    end
        
    %visualize what I am doing
    figure(200)
    plot(timeDay,szCorr*100,'r.','MarkerSize',20)
    hold on;
    grid on;
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1);
    set(gcf,'position',get(0,'screensize'));
    vline(timeDay)
    vline(szType.szT.test(1)/60/60/24,'b')
    hline(60,'m');
    ylabel('Correlation (%)')
    xlabel('Time (days)')
    title('Correlation of Predicted Time to Seizure for Each Seizure')   
%     plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\szCorrPlots_9-9\' pt{i} 'corrBySzretrained'];
    plotName = ['H:\jaredwil\Lasso Results\retrainStuff\' pt{i} 'corrBySzretrained'];
    savefig(plotName)
    saveas(gcf,plotName,'jpg')
%     close all;
end