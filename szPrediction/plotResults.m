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
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\5minFeats'))  %this is where .mat file are contained on local comp
% addpath(genpath('H:\jaredwil\szPred_feats\Win_5min'))  %this is where .mat file are contained on local comp

set(0,'DefaultTextInterpreter','none');
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

%%
%begin function
%loop through all pts
for i = 1:numel(pt);  %%%%%TEMPORARY for debug
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

    featLab = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Res LassoLib\' pt{i} '_bestLasso5.mat'];
%     featLab = ['H:\jaredwil\Lasso Results\Res 8-29\' pt{i} '_bestLasso5.mat'];

    load(featLab);


    f = lassoRes.coef.totmin;
    tInt = lassoRes.int(3);
    tInt2 = repmat(tInt, size(allFeat,1),1);


    predLabOrg = allFeat*f + tInt2;

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
        szCorr(sz) = corr(allLabs(uIdx == sz),predLabOrg(uIdx == sz));
        %random correlation

        
        timeDay(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(allSzLab(uIdx == sz,1));
        
        trainSz(sz) = sum((timeSec(sz) == szType.szT.train(:,1)));
    end
    
    %eliminate sz that do not have many observationsbecause
    %it messes up the null hypothesis
    idxElim = (winDur_sz < 60);
    szCorrHist = szCorr;
    szCorrHist(idxElim)        = [];
    timeDay(idxElim)       = [];
%     timeSec(idxElim)       = [];
    trainSz(idxElim)       = [];
    
    
% Original corr by sz has been found now threshold and keep sz's with high
% correlation and plot results
    corrTh = .65; %Null Hyptho Test

    corrLab = zeros(size(allLabs));
    isTest = zeros(size(allLabs));
    retrainSz = zeros(numel(numSz),1);

    trainSztimes = szType.szT.train(:,1);
    trainSztimes = trainSztimes(trainSz == 1);
    
    testSztimes = szType.szT.test(:,1);
    % trainSztimes = trainSztimes(1:round(length(trainSztimes) * (0.4/0.7)));  %This reduces the training set from 70% to 30%

    for sz = 1:numel(numSz);
        corrLab(uIdx == sz) = szCorr(sz);                                   %repeat correlation for a label
        isTest(uIdx == sz) = (sum((timeSec(sz) == testSztimes)));  %sz was orignally trained on
        retrainSz(sz) = sum((timeSec(sz) == trainSztimes)) & (szCorr(sz) > corrTh); %the retrain sz index is based on sz corr and time 
    end

    %these are the feature idx to train a new model on and it is important to
    %note this because they both are inside the new "training set" and are also
    %withing the correlation threshold for retraining
    idxplotTest = (corrLab > corrTh) & isTest; 

    plotFeats  = allFeat(idxplotTest,:);
    plotTestLabels = allLabs(idxplotTest);
    plotPred = predLabOrg(idxplotTest);
    
    %visual assesment
    time5 = linspace(5,length(plotPred)*5,length(plotPred))/60;
    figure(1)
    set(gca,'FontSize',15);
    set(gca,'LineWidth',2);
    set(gcf,'Position',get(0,'Screensize')); 
    %make background white
    set(gcf,'Color','w');
    plot(time5, plotTestLabels/60/60);
    hold on;
    plot(time5, smooth(plotPred/60/60,5));
%     title(['Time to Seizure Visualization of Predictions with Better Than Random Corr.: pt. ' pt{i}])
    title(['Time to Seizure Visualization of Predictions with Corr. > 65%: pt. ' pt{i}])

    xlabel('Time (Hours)')
    ylabel('Time to Seizure (Hours)')
    legend('Test','Prediction','location','best')
    plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\Meeting Results 9-20\' pt{i} 'highCorr_Results'];
    saveas(gcf,plotName,'jpg')
    savefig(gcf,plotName)
    
end