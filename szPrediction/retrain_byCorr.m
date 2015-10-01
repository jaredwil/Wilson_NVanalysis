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
[p, poolsize ] = initParPool( 10 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lassoTime = zeros(1,numel(pt));
%%
%begin function
%loop through all pts
for i = 6%1:numel(pt);  %%%%%TEMPORARY for debug
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
    
    
%%%%%%%%%%%%%% REMOVE SZ THAT ARE NOT CCS (CLINICALLY CONFIRMERD) %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ONLY LOOK AT CCS
    allFeat(allSzLab(:,2)~=1,:) = [];
    allLabs(allSzLab(:,2)~=1) = [];
    allSzLab(allSzLab(:,2)~=1,:) = [];

    %load feature weights
    % load('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\resultsV1\NVC1001_23_005_bestLasso_CCS.mat')

    % featLab = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Res LassoLib\' pt{i} '_bestLasso5.mat'];
    featLab = ['H:\jaredwil\Lasso Results\Res 8-29\' pt{i} '_bestLasso5.mat'];

    load(featLab);


    f = lassoRes.coef.totmin;
    tInt = lassoRes.int(3);
    tInt2 = repmat(tInt, size(allFeat,1),1);


    predLabOrg = allFeat*f + tInt2;

    %visual assesment
    figure(1)
    plot(allLabs);
    hold on;
    plot(smooth(predLabOrg,5));


    %Create a random label matrix mx1000 where m is the number of labels
    randLab = rand(length(predLabOrg),1000)*7200;
    permLab = rand(length(predLabOrg),1000)*7200;


    %average feature across sz do pca then cluster. 
    [numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 
    
    szCorr = zeros(numel(numSz),1);
    randCorr = zeros(numel(numSz),1000);  %correlation of random lables
    permCorr = zeros(numel(numSz),1000);  %correlation of permuted pred lables

    timeDay   = zeros(numel(numSz),1);
    timeSec   = zeros(numel(numSz),1);
    trainSz   = zeros(numel(numSz),1);
    winDur_sz = zeros(numel(numSz),1);
    for sz = 1:numel(numSz);
        
        winDur_sz(sz) = length(allLabs(uIdx == sz));
        
        %create permutation labels rather than random labels
        for permIt = 1:1000
            tmp = predLabOrg(uIdx == sz);
            permord = randperm(winDur_sz(sz));
            permLab(uIdx == sz,permIt) = tmp(permord);
        end
        
%       szAvg_feat(sz,:) = [mean(allFeat(uIdx == sz,:),1) szTimes(sz,1)];
        szCorr(sz) = corr(allLabs(uIdx == sz),predLabOrg(uIdx == sz));
        %random correlation
        randCorr(sz,:) = corr(allLabs(uIdx == sz),randLab(uIdx == sz,:));
        %permuted predicted labels for each sz rather than random
        %generation
        permCorr(sz,:) = corr(allLabs(uIdx == sz),permLab(uIdx == sz,:));
        
        timeDay(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(allSzLab(uIdx == sz,1));
        
        trainSz(sz) = sum((timeSec(sz) == szType.szT.train(:,1)));
    end
    
    %eliminate sz that do not have many observationsbecause
    %it messes up the null hypothesis
    % this does two things not really sure if I should be doing both at the
    % same time The first part eliminates short sz to remove false high sz
    % correlations, the second eliminates sz in the first 15 days. Features
    % Normalization
    idxElim = (winDur_sz < 6) | (timeDay < 15);  
    szCorrHist = szCorr;
    szCorrHist(idxElim)        = [];
    randCorr(idxElim,:)    = [];
    permCorr(idxElim,:)    = [];

    timeDay(idxElim)       = [];
%     timeSec(idxElim)       = [];
    trainSz(idxElim)       = [];
    
%     rCorV = randCorr(:);
    rCorV = permCorr(:);

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
    vline(szType.szT.test(1)/60/60/24,'b','Tr/Tst Divide')
    hline(corr5*100,'m');
    ylabel('Correlation (%)')
    xlabel('Time (days)')
    title('Correlation of Predicted Time to Seizure for Each Seizure')
    plotName = ['H:\jaredwil\Lasso Results\retrainStuff\' pt{i} 'corrBySzLib'];
    saveas(gcf,plotName,'jpg')
    
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
%     trainSztimes = trainSztimes(trainSz == 1);
%     trainSztimes(trainSz == 1) = 0;

    testSztimes = szType.szT.test(:,1);
    % trainSztimes = trainSztimes(1:round(length(trainSztimes) * (0.4/0.7)));  %This reduces the training set from 70% to 30%

    for sz = 1:numel(numSz);
        corrLab(uIdx == sz) = szCorr(sz);                                   %repeat correlation for a label
        isTrain(uIdx == sz) = sum((timeSec(sz) == trainSztimes));  %sz was orignally trained on
        retrainSz(sz) = any((timeSec(sz) == trainSztimes)) & (szCorr(sz) > corrTh); %the retrain sz index is based on sz corr and time 
    end
    
    retrainSz(idxElim) = [];
    figure(2000)
    plot(timeDay(retrainSz == 1),szCorrHist(retrainSz == 1)*100,'b.','MarkerSize',20)  %Test sz
    hold on;
    grid on;
    plot(timeDay(retrainSz ~= 1),szCorrHist(retrainSz ~= 1)*100,'r.','MarkerSize',20)  %Test sz
    title('Retraining Distribution')
    vline(15, 'm', 'Day 15')
    vline(timeDay)
    vline(szType.szT.test(1)/60/60/24,'b','tr/tst')
    hline(corr5*100,'m');
    legend('Retrain','Not Retrain')

    %these are the feature idx to train a new model on and it is important to
    %note this because they both are inside the new "training set" and are also
    %withing the correlation threshold for retraining
    idxTrain = (corrLab > corrTh) & isTrain; 

    retrainFeats  = allFeat(idxTrain,:);
    retrainLabels = allLabs(idxTrain);


    if(isempty(retrainFeats))
       disp(['pt ' pt{i} ' does not have enough data for retrain'])
       continue; 
    end


    numLam = 75;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method Using Lasso Lib
    disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
    lambda   = linspace(1e4,1e6,numLam);
%   lambda   = logspace(4,6,numLam);
    numFolds = 10;
    numFolds = 3;
    modCh = 'se';
    try
        [ f, tInt, numFeats, solnLambda, lassoTime(i) ] = trainLassoLib( retrainFeats,retrainLabels,numFolds,lambda, modCh, pt{i});
    catch
        disp('Could Not Find Solution')
    end
    
    disp(['DONE Training Lasso Model on Patient: ' pt{i}])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % %Other method -- Uses built in Lasso Soln. Function
    % numFolds = 3;
    % modCh = 'min';
    % [ f, tInt, numFeats, lassoTime(i) ] = trainLassoReg( retrainFeats,retrainLabels,numFolds,lambda, modCh );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %NOW REDO THE FIRST PART WITH THE NEW STUFFS

    % 
    % f = bestLasso_test;
    % tInt = int_test;
    % f = bestLasso_corr;
    % tInt = bestInt_corr;
    tInt2 = repmat(tInt, size(allFeat,1),1);


    predRe = allFeat*f + tInt2;

    %visual assesment
    figure(100)
    plot(allLabs);
    hold on;
    plot(smooth(predLabOrg,5));
    plot(smooth(predRe,5));
    xlabel('Time (kinda)')
    ylabel('Time to Seizure (seconds)')
    legend('True','Original','Retrained')
    title(['Compare Resulting Labels for Paitient ' pt{i}]);
    plotName = ['H:\jaredwil\Lasso Results\retrainStuff\' pt{i} 'compareResults'];
    saveas(gcf,plotName,'jpg')
    
    
    %average feature across sz do pca then cluster. 
    [numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 
    
    szCorr = zeros(numel(numSz),1);
    timeDay = zeros(numel(numSz),1);
    timeSec = zeros(numel(numSz),1);
    for sz = 1:numel(numSz);
%         szAvg_feat(sz,:) = [mean(allFeat(uIdx == sz,:),1) szTimes(sz,1)];
        szCorr(sz) = corr(allLabs(uIdx == sz),predRe(uIdx == sz));
        timeDay(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(allSzLab(uIdx == sz,1));
    end
    %eliminate the same sz as before
%     szCorrHist = szCorr;
    szCorr(idxElim)        = [];
    timeDay(idxElim)       = [];
    
    %save a struct to look at data later
    retrainInfo = struct('pred', predRe, 'allLabels',allLabs,'szCorr',szCorr, 'f',f,'tInt',tInt,'lambda',solnLambda);
    infoLabel = ['H:\jaredwil\Lasso Results\retrainStuff\' pt{i} 'retrainInfoLib.mat'];
    save(infoLabel,'retrainInfo')
    
    %RETRAINED HIST
    %make a histogram of the resulting test sz to see how many are better
    %than random with confidence (p < 0.05).
    %only make a histogram of the sz not retrained on & in test sz range
    testSzCorr = szCorr(retrainSz ~= 1 & (timeDay > szType.szT.train(end,1)/60/60/24));
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
    saveas(gcf,plotName,'jpg')
end