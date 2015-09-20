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
% [p, poolsize ] = initParPool( 10 );

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
for i = 5:numel(pt);  %%%%%TEMPORARY for debug
     close all;

%   label = [pt{i} '_szPred_30secFeats.mat'];
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


    %Create a random label matrix mx1000 where m is the number of labels
    randLab = rand(length(predLab),1000)*7200;
    

    %average feature across sz do pca then cluster. 
    [numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 
    
    szCorr = zeros(numel(numSz),1);
    randCorr = zeros(numel(numSz),1000);

    timeDay   = zeros(numel(numSz),1);
    timeSec   = zeros(numel(numSz),1);
    trainSz   = zeros(numel(numSz),1);
    winDur_sz = zeros(numel(numSz),1);
    for sz = 1:numel(numSz);
        
        winDur_sz(sz) = length(allLabs(uIdx == sz));
%       szAvg_feat(sz,:) = [mean(allFeat(uIdx == sz,:),1) szTimes(sz,1)];
        szCorr(sz) = corr(allLabs(uIdx == sz),predLab(uIdx == sz));
        %random correlation
        randCorr(sz,:) = corr(allLabs(uIdx == sz),randLab(uIdx == sz,:));

        
        timeDay(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(allSzLab(uIdx == sz,1));
        
        trainSz(sz) = sum((timeSec(sz) == szType.szT.train(:,1)));
    end
    
    %eliminate sz that do not have many observationsbecause
    %it messes up the null hypothesis
    idxElim = (winDur_sz < 6);
    szCorrHist = szCorr;
    szCorrHist(idxElim)        = [];
    randCorr(idxElim,:)    = [];
    timeDay(idxElim)       = [];
%     timeSec(idxElim)       = [];
    trainSz(idxElim)       = [];
    
    rCorV = randCorr(:);
%find corr where p < 0.05
    tmp = .3:0.0001:.5;
    pVals = [];
    for pIdx = 1:length(tmp)
        pVals(pIdx) = (sum(rCorV > tmp(pIdx))/size(rCorV,1))*100 ;
    end
    findp = abs(pVals - 5);
    [idx idx] = min(findp);
    corr5 = tmp(idx)
    
    figure(22)
    randCorr = mean(randCorr,2);
    hist(rCorV,1000);
    hold on;
%     corr5 = 0.35517; 
    p5 = (sum(rCorV > corr5)/size(rCorV,1))*100 
    vline(corr5,'r:',['p = 0.5, Null Corr = ' num2str(corr5)])
    
    %This was found using the below script
    tmp = .3:0.0001:.4;
    
    pVals = [];
    for pIdx = 1:length(tmp)
        pVals(pIdx) = (sum(rCorV > tmp(pIdx))/size(rCorV,1))*100 ;
    end

    
    %visualize what I am doing
    figure(2)
    plot(timeDay(trainSz == 1),szCorrHist(trainSz == 1)*100,'b.','MarkerSize',20)  %Test sz
    hold on;
    grid on;
    plot(timeDay(trainSz ~= 1),szCorrHist(trainSz ~= 1)*100,'r.','MarkerSize',20)  %Test sz
%     plot(timeDay,randCorr*100,'b*','MarkerSize',20)  %this doesn't really
%     tell me anything
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1);
    set(gcf,'position',get(0,'screensize'));
    vline(timeDay)
    vline(szType.szT.test(1)/60/60/24,'b')
    hline(corr5*100,'m');
    ylabel('Correlation (%)')
    xlabel('Time (days)')
    title('Correlation of Predicted Time to Seizure for Each Seizure')
    
    %only interested in test sz correlation
    testSzCorr = szCorrHist(trainSz ~= 1);
%Original HIST
%make a histogram of the resulting test sz to see how many are better
%than random with confidence (p < 0.05).
    %first plot hist or test corr
    histRes = 20;
    histBins = linspace(-1,1,histRes);
    figure(3)
    [C1, C1x] = hist(testSzCorr,histBins);
    C1 = C1./sum(C1);
    bar(C1x,C1)
    hold on;
    %plot null distribution for comparison
    C2x = linspace(-1,1,histRes);
    C2 = hist(rCorV,histRes);
    C2 = C2./sum(C2);
    plot(C2x,C2)
    vline(corr5,'r:','p = 0.05')
    xlabel('Test Correlation')
    ylabel('Density (%)')
    title('Distribution of Testing Correlations Compared to Null Hypothesis (Original)')
    legend('Test','Null')
    xlim([-1 1])
    histName = ['H:\jaredwil\Lasso Results\retrainStuff\' pt{i} 'histOriginal'];
    savefig(histName)
    saveas(gcf,histName,'jpg')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Original corr by sz has been found now threshold and keep sz's with high
% correlation then retrain

corrTh = corr5; %Null Hyptho Test

corrLab = zeros(size(allLabs));
isTrain = zeros(size(allLabs));
retrainSz = zeros(numel(numSz),1);

trainSztimes = szType.szT.train(:,1);
trainSztimes = trainSztimes(trainSz == 1);
% trainSztimes = trainSztimes(1:round(length(trainSztimes) * (0.4/0.7)));  %This reduces the training set from 70% to 30%

for sz = 1:numel(numSz);
    corrLab(uIdx == sz) = szCorr(sz);                                   %repeat correlation for a label
    isTrain(uIdx == sz) = sum((timeSec(sz) == trainSztimes));  %sz was orignally trained on
    retrainSz(sz) = sum((timeSec(sz) == trainSztimes)) & (szCorr(sz) > corrTh); %the retrain sz index is based on sz corr and time 
end

%these are the feature idx to train a new model on and it is important to
%note this because they both are inside the new "training set" and are also
%withing the correlation threshold for retraining
idxTrain = (corrLab > corrTh) & isTrain; 

retrainFeats  = allFeat(idxTrain,:);
retrainLabels = allLabs(idxTrain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method Using Lasso Lib
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
plot(lambda, tLasso_MSE,'r.-')

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
h = plot(dzT(corrIdx),tLasso_coor(corrIdx)*100,'bs','LineWidth',4);
% h = vline(dzT(seIdx),'b:');
plotName = ['H:\jaredwil\Lasso Results\' pt{i} '_retrainCorrLib'];
saveas(h,plotName,'jpg')

testIdx = minIdx - 4;

[bestLasso_test, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(testIdx),'maxIter',50000);
[bestLasso_corr, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(corrIdx),'maxIter',50000);
[bestLasso_Min, ~, ~] = LassoBlockCoordinate(retrainFeats,retrainLabels,lambda(minIdx),'maxIter',50000);


numFeats_test = sum(bestLasso_test ~= 0);
numFeats_corr = sum(bestLasso_corr ~= 0);
numFeats_Min = sum(bestLasso_Min ~= 0);

bestInt_corr = mean(retrainLabels) -  mean(retrainFeats*bestLasso_test);
bestInt_test = mean(retrainLabels) -  mean(retrainFeats*bestLasso_corr);
bestInt_Min  = mean(retrainLabels) -  mean(retrainFeats*bestLasso_Min);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %Other method
%split up for a simple 2 fold CV
% tr = length(retrainLabels)*0.7;
% evalFeats = retrainFeats(tr+1:end,:);
% evalLabels = retrainLabels(tr+1:end);
% retrainFeats = retrainFeats(1:tr,:);
% retrainLabels = retrainLabels(1:tr); 
% 
% tic;
% disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
% opts = statset('UseParallel',true);  %Do dis shit in parallel!!!  use number of availble workers for CV
% [fLasso, tINFO] = lasso(retrainFeats, retrainLabels,'CV',5,'Options',opts);
% lassoTime(i) = toc;
% %BOOM!!!
% disp(['DONE Training Lasso Model on Patient: ' pt{i}])
% dzT = tINFO.DF; 
% tInt = tINFO.Intercept;
% % xAlpha = xINFO.Alpha;
% tInt2 = repmat(tInt, size(retrainFeats,1),1);
% utLasso = retrainFeats*fLasso + tInt2;
% tTrain = repmat(retrainLabels,1,100);
% tLasso_coor = corr(tTrain,utLasso);
% tLasso_coor = tLasso_coor(1,:);
% 
% figure(11)
% plot(dzT,tLasso_coor*100,'r.');
% xlabel('Number of Non-Zero Feature Weights in Lasso Model')
% ylabel('Resulting Test Correlation (%)')
% title('TRAIN COULD BE MISLEADING')
% grid on;
% hold on;
% 
% tInt2 = repmat(tInt, size(evalFeats,1),1);
% utLasso = evalFeats*fLasso + tInt2;
% tTrain = repmat(evalLabels,1,100);
% tLasso_coor = corr(tTrain,utLasso);
% tLasso_coor = tLasso_coor(1,:);
% 
% minIdx = tINFO.IndexMinMSE;
% seIdx  = tINFO.Index1SE;
% corrIdx = tLasso_coor == max(tLasso_coor);
% 
% figure(12)
% plot(dzT,tLasso_coor*100,'r.');
% xlabel('Number of Non-Zero Feature Weights in Lasso Model')
% ylabel('Resulting Test Correlation (%)')
% title('EVALUATION')
% grid on;
% hold on;
% h = plot(dzT(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
% h = vline(dzT(minIdx),'g:');
% h = plot(dzT(seIdx),tLasso_coor(seIdx)*100,'bs','LineWidth',4);
% h = vline(dzT(seIdx),'b:');
% 
% figure(14)
% h2 = lassoPlot(fLasso,tINFO,'PlotType','CV');
% % Use a log scale for MSE to see small MSE values better
% % set(gca,'YScale','log');
% % plotName = ['H:\jaredwil\Lasso Results\matLasso_history\' pt{i} '_mseRes_CCS'];
% 
% bestLasso_corr = fLasso(:,corrIdx);
% bestLasso_1SE = fLasso(:,seIdx);
% bestLasso_Min = fLasso(:,minIdx);
% 
% bestInt_corr = tInt(corrIdx); 
% bestInt_1SE = tInt(seIdx); 
% bestInt_Min = tInt(minIdx); 
% 
% numFeats_corr = dzT(corrIdx); 
% numFeats_1SE = dzT(seIdx); 
% numFeats_Min = dzT(minIdx);
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
    
    %save a struct to look at data later
    retrainInfo = struct('pred', predLab, 'allLabels',allLabs,'szCorr',szCorr, 'f',f,'tInt',tInt);
    infoLabel = ['H:\jaredwil\Lasso Results\retrainStuff\' pt{i} 'retrainInfoLib.mat'];
    save(infoLabel,'retrainInfo')
    
    %RETRAINED HIST
    %make a histogram of the resulting test sz to see how many are better
    %than random with confidence (p < 0.05).
    testSzCorr = szCorr(retrainSz ~= 1);
    histRes = 20;
    histBins = linspace(-1,1,histRes);
    figure(199)
    [C1, C1x] = hist(testSzCorr,histBins);
    C1 = C1./sum(C1);
    bar(C1x,C1)
    hold on;
    C2x = linspace(-1,1,histRes);
    C2 = hist(rCorV,histRes);
    C2 = C2./sum(C2);
    plot(C2x,C2)
    vline(corr5,'r:','p = 0.05')
    xlabel('Test Correlation')
    ylabel('Density (%)')
    title('Distribution of Testing Correlations Compared to Null Hypothesis (RETRAINED)')
    legend('Test','Null')
    histName = ['H:\jaredwil\Lasso Results\retrainStuff\' pt{i} 'histRetrainedLib'];
    savefig(histName)
    saveas(gcf,histName,'jpg')
    
    
    %visualize what I am doing
    figure(200)
    plot(timeDay(retrainSz == 1),szCorr(retrainSz == 1)*100,'b.','MarkerSize',20)
    hold on;
    grid on;
    plot(timeDay(retrainSz ~= 1),szCorr(retrainSz ~= 1)*100,'r.','MarkerSize',20)
    set(gca,'FontSize',12);
    set(gca,'LineWidth',1);
    set(gcf,'position',get(0,'screensize'));
    vline(timeDay)
    vline(trainSztimes(end)/60/60/24 + 1,'b','tr/test divide')
    hline(corr5*100,'m','Null Dist.');
    ylabel('Correlation (%)')
    xlabel('Time (days)')
    title('Correlation of Predicted Time to Seizure for Each Seizure')   
    legend('Train','Test')
%     plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\szCorrPlots_9-9\' pt{i} 'corrBySzretrained'];
    plotName = ['H:\jaredwil\Lasso Results\retrainStuff\' pt{i} 'corrBySzretrainedLib'];
    savefig(plotName)
    saveas(gcf,plotName,'jpg')
end