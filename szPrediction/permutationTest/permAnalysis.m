%Permutation Test Analysis
%Creat permutation of labels and train a new model to see if original
%results are signif. (p < 0.05)

%%
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp



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
idxSuccess = cell(numel(pt),1);
numSuc     = zeros(numel(pt),2);
perSuc     = zeros(numel(pt),1);
for i = 1:numel(pt)
    
    try
        load(['H:\jaredwil\Lasso Results\permutationTest\' pt{i} '\permInfo.mat'])
    catch
        disp(['NO INFO FOR PT. ' pt{i}])
        continue;
    end
   
    close all;
     clc;
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
    
    trainFeats =  bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);
    testFeats =  bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);
        
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

    fCorr = lassoRes.coef.corr;
    tIntCorr = lassoRes.int(1);
    tIntCorr2 = repmat(tIntCorr, size(allFeat,1),1);
    predCorr = allFeat*fCorr + tIntCorr2;

    
    
    fMin = lassoRes.coef.min;
    tIntMin = lassoRes.int(2);
    tIntMin2 = repmat(tIntMin, size(allFeat,1),1);
    predMin = allFeat*fMin + tIntMin2;

    f = lassoRes.coef.totmin;
    tInt = lassoRes.int(3);
    tInt2 = repmat(tInt, size(allFeat,1),1);
    testInt = repmat(tInt, size(testFeats,1),1);

    pred = allFeat*f + tInt2;
    predTest = testFeats*f + testInt;


    %average feature across sz do pca then cluster. 
    [numSz_train,~, uIdx] = unique(trainSzT(:,1),'stable'); 
    [numSz_test,~, uIdxTest] = unique(testSzT(:,1),'stable'); 
    
    timeDay   = zeros(numel(numSz_train),1);
    timeSec   = zeros(numel(numSz_train),1);
    winDur_sz = zeros(numel(numSz_train),1);
    for sz = 1:numel(numSz_train);
        
        winDur_sz(sz) = length(trainLabels(uIdx == sz));
                
        timeDay(sz) = mode(trainSzT(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(trainSzT(uIdx == sz,1));
        
    end
        
%%%%%%compute test correlation to compare against null model


    testInt = repmat(tInt, size(testFeats,1),1);
    predTest = testFeats*f + testInt;
    testCorr = zeros(numel(numSz_test,1));
    winDur_test = zeros(numel(numSz_test),1);
    testSz_type = zeros(numel(numSz_test),1);
    for sz = 1:numel(numSz_test);
        winDur_test(sz) = length(testLabels(uIdxTest == sz));
        testSz_type(sz) = mode(testSzT(uIdxTest == sz,2));

        testCorr(sz,:) = corr(predTest(uIdxTest == sz),testLabels(uIdxTest == sz));
    end
    idxElim = (winDur_test < 6) | (testSz_type ~= 1); %remove if short Sz or not CCS sz.
    testCorr(idxElim) = [];
   
%%%%%%%%%%%%%%%%% This is where the magic happens %%%%%%%%%%%%%%%%%%%%%%%%%
%% check if the test corr was higher than the null dist.
    tmpSuc = zeros(length(testCorr),1);
    for sz = 1:length(testCorr);
        tmpSuc(sz) = testCorr(sz) >= permInfo(sz,3).nullCorr;
    end
    idxSuccess{i} = tmpSuc;
    
    numSuc(i,1) = sum(tmpSuc);
    numSuc(i,2) = length(testCorr);
    if(isempty(tmpSuc) == 0)
        perSuc(i) = sum(tmpSuc)/length(tmpSuc);
    end
end
permRes = struct('precentValid',[],'numberValid',[],'validIdx',[],'overallValid',sum(numSuc(:,1))/sum(numSuc(:,2)));
permRes.precentValid = perSuc;
permRes.numberValid = numSuc;
permRes.validIdx = idxSuccess;

% permCell = {perSuc,numSuc,idxSuccess,sum(numSuc(:,1))/sum(numSuc(:,2))};
% permRes = struct('precentValid',perSuc,'numberValid',numSuc,'validIdx',idxSuccess,'overallValid',sum(numSuc(:,1))/sum(numSuc(:,2)));
infoLabel = 'H:\jaredwil\Lasso Results\permutationTest\permResults.mat';
save(infoLabel,'permRes')

M = [numSuc(:,2) numSuc(:,1) perSuc];
csvwrite('H:\jaredwil\Lasso Results\permutationTest\permResults.csv',M)
    
    
