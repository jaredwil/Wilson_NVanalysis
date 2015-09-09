%Jared Wilson
%Look at correlation by seizure


%%
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\30secFeats'))  %this is where .mat file are contained on local comp
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\5minFeats'))  %this is where .mat file are contained on local comp


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
for i = 12;  %%%%%TEMPORARY for debug
       
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

%%%%%%%%%%%%%%%%%%%%%%%RESHAPE FEATURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1:2; %number of history samples
samples = size(allFeat,1);
features = size(allFeat,2);
startSAMP = max(N)+1;
M = (samples-length(N));  %number of time bins
%create 'R' matrix for linear regression algorithm
r = zeros(M, features*length(N)+1);
for resIdx = 1:M
    temp = allFeat(startSAMP + (resIdx-1) - N,:);   %temp is a temporary matrix    
    r(resIdx,:) = [1 temp(:)'];
end

allLabs = allLabs(startSAMP:end,:);
allSzLab = allSzLab(startSAMP:end,:);

[~,IdxRemove] = findpeaks(allLabs);

r(IdxRemove,:) = [];
allSzLab(IdxRemove,:) = [];
allLabs(IdxRemove,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%load feature weights
load('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\resultsV1\NVC1001_23_005_bestLasso_CCS.mat')

f = lassoRes.coef.min;
tInt = lassoRes.int(2);
tInt2 = repmat(tInt, size(r,1),1);


predLab = r*f + tInt2;

%visual assesment
figure(1)
plot(allLabs);
hold on;
plot(smooth(predLab,5));


    %average feature across sz do pca then cluster. 
    [numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 
    
    szCorr = zeros(numel(numSz),1);
    time = zeros(numel(numSz),1);

    for sz = 1:numel(numSz);
%         szAvg_feat(sz,:) = [mean(allFeat(uIdx == sz,:),1) szTimes(sz,1)];
        szCorr(sz) = corr(allLabs(uIdx == sz),predLab(uIdx == sz));
        time(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
    end
    
    figure(2)
    plot(time,szCorr,'r.','MarkerSize',20)
    hold on;
    grid on;
    vline(time)

end