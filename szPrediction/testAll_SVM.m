%Sz Prediction Pipeline v1 -- use of functions as much as possible
%Jared Wilson
%7/24/2015

% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))

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
% i = 12;  %%%%%TEMPORARY for debug
allFeats_train  = cell(numel(pt),1);
allLab_train    = cell(numel(pt),1);
model           = cell(numel(pt),1);
allFeats_test   = cell(numel(pt),1);
allLab_test     = cell(numel(pt),1);
accuracy        = cell(numel(pt),1);
predLab         = cell(numel(pt),1);


%create log file
fileID = fopen('svmLog.txt','wt');
fclose(fileID);
for i = 1:numel(pt)

    %Get train features and training labels (lables -> minutes to sz)
    %features extracted defined in szPred_winFeatExt <- change this file
    [trainFeats, trainLabels ] = szPred_train(pt{i}, usernm, pswdBin, trPct, winLen, winDisp, szHorizon);

    %Skip if there are no train features meaning no sz annotation
    if(isempty(trainFeats))
        disp([pt{i} ' SKIPPPED'])
       continue; 
    end

    %normalize features
    avgFeats = mean(trainFeats,1);
    stdFeats = std(trainFeats,[],1);
    trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);

    %SVM Labels
    svmLabels = trainLabels;
    svmLabels(svmLabels >= 0 & svmLabels < szHorizon*60*60) = 1; %preictal
    svmLabels(svmLabels > szHorizon*60*60) = 0; %interictal

    %keep everything so we can go back and check it
    allFeats_train{i} = trainFeats;
    allLab_train{i} = trainLabels;

    %%
    %HERE IS WHERE YOU WOULD TRAIN YOUR MODEL
    % Train SVM model
    model{i} = svmtrain(svmLabels, trainFeats, '-t 0 -b 1 -c 1'); 

    %%
    % Get TEST FEATS
    %get the training features and labels
    [testFeats, testLabels] = szPred_test(pt{i}, usernm, pswdBin, trPct, winLen, winDisp, szHorizon);
    %normalize test feats
    testFeats = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);


    %%
    %make predictions based on training model

    %SVM
    svmTestLabels = testLabels;
    svmTestLabels(svmTestLabels >= 0 & svmTestLabels < szHorizon*60*60) = 1; %preictal
    svmTestLabels(svmTestLabels > szHorizon*60*60) = 0; %interictal

    [predLab{i}, accuracy{i}, prob_estimates] = svmpredict(svmTestLabels, testFeats, model{i}, '-b 1');

    allFeats_test{i} = testFeats;
    allLab_test{i} = testLabels;
    
    %write a log file
    fileID = fopen('svmLog.txt','w');
    nbytes = fprintf(fileID,'%f\n',i);
    fclose(fileID);
    
end