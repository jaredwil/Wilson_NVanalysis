%Jared Wilson
%8/30/2015

%clust seizure using k-means

%%
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp

%%
% Define algorithm specifics 
usernm = 'jaredwil';
pswdBin = 'jar_ieeglogin.bin';
trPct = 0.7;
winLen = 30;
winDisp = 30;
szHorizon = 1; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%%
%begin function
%loop through all pts
for i = 12 ;  %%%%%TEMPORARY for debug
       
    label = [pt{i} '_szPred_30secFeats.mat'];
    tylabel = [pt{i} '_szTypeLabels.mat'];

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
    
    
    allFeat = [trainFeats; testFeats];
    allSzLab  = [trainSzT;testSzT];
    
    
    %szTimeline
    figure(1)
    vline(szTimes(:,1)/60/60/24)
    
    %average feature across sz do pca then cluster. 
    [~,~, uIdx] = unique(allSzLab(:,1),'stable');                    
    
    
    
    
    
    
    
    




end