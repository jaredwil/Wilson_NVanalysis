%Script to test grabbing interictal data randomly vs consecutivly
%Jared Wilson
%7/21/2015

% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))

pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

i = 12;  %%%%%TEMPORARY

session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
fs = session.data.sampleRate;               %Find sampling Rate
 
%get ALL seziure start and end times from anotations
[startT, endT] = getSzTime(session);

%split up seziures into test and training sets
%precent of total seizures to train on;
trPct = 0.4;

numTr = floor(length(startT)*trPct);  %number of training sz
numTest = length(startT) - numTr;     %number of testing sz

trainST = startT(1:numTr);
trainET = endT(1:numTr);
trInterIct_endT = trainET(end);   %use the final sz end time of train set 
                                  %as the end training interictal search
testST = startT(numTr+1:end);
testET = endT(numTr+1:end);

szHorizon = 2; %in hours
winLen = 30;
winDisp = 30;
%Get training Data
%Sz horizon data
[train, predIdx] = getWinFeats_Train(session, trainST, trainET, szHorizon, winLen, winDisp);
%interictal data retrieved from the first available blocks
[trainInt] = Interict_train( session, predIdx, winLen, winDisp);

%isolate time lables and training feats
szhorzTrain = train(:,2:end);
szhorzLabels = train(:,1);
intLables = ones(size(trainInt,1),1)*(szHorizon*60*60 + winLen);

trainLabels = [szhorzLabels; intLables];
trainFeats = [szhorzTrain;trainInt];

%normalize features
avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);
trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);

%%%%GET interictal training data from random blocks rather than the first
%%%%few hours
[trainRandInt] = InterictRand_train( session, predIdx, winLen, winDisp);
intRandLables = ones(size(trainRandInt,1),1)*(szHorizon*60*60 + winLen);

trainLabels_rand = [szhorzLabels; intRandLables];
trainFeats_rand = [szhorzTrain;trainRandInt];

%normalize features
avgFeats_rand = mean(trainFeats_rand,1);
stdFeats_rand = std(trainFeats_rand,[],1);
trainFeats_rand = bsxfun(@rdivide, bsxfun(@minus,trainFeats_rand,avgFeats_rand), stdFeats_rand);

%%
%%%%%%%%%%%%%%%%%%%%% TRAINING SOME MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try this logistic regresion PDF thing
% lrLables = categorical(trainLabels/30); 
% 
% B = mnrfit(trainFeats,lrLables);

%for K-nn
%try a very simple classification classifying 2hr (2) vs 1hr (1) vs
%interictal (3)
knnLabels = trainLabels;
knnLabels(knnLabels >= 0 & knnLabels < 1*60*60) = 1; %between 2-1 hours
knnLabels(knnLabels >= 1*60*60 & knnLabels < szHorizon*60*60) = 2; %between 2-1 hours
knnLabels(knnLabels > szHorizon*60*60) = 3; %interictal

knnLabels_rand = trainLabels_rand;
knnLabels_rand(knnLabels_rand >= 0 & knnLabels_rand < 1*60*60) = 1; %between 2-1 hours
knnLabels_rand(knnLabels_rand >= 1*60*60 & knnLabels_rand < szHorizon*60*60) = 2; %between 2-1 hours
knnLabels_rand(knnLabels_rand > szHorizon*60*60) = 3; %interictal

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get seizure testing data
[test, testpredIdx] = getSzFeats(session, testST, testET, szHorizon, winLen, winDisp);
%get interictal testing data
[testInt] = Interict_test(session, testpredIdx, winLen, winDisp);

%create a single test matrix
szhorzTest = test(:,2:end);
szhorzTsLabels = test(:,1);
testIntLables = ones(size(testInt,1),1)*(szHorizon*60*60 + winLen);

testLabels = [szhorzTsLabels; testIntLables];
testFeats = [szhorzTest;testInt];

%normalize features based on training mean and std
testFeats_nr = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);

testFeats_rand = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats_rand), stdFeats_rand);

%CREATE LABELS  classification classifying 1hr (1) vs 2hr (2) vs interictal (3)
testLabels(testLabels >= 0 & testLabels < 1*60*60) = 1; %between 2-1 hours
testLabels(testLabels >= 1*60*60 & testLabels < szHorizon*60*60) = 2; %between 2-1 hours
testLabels(testLabels > szHorizon*60*60) = 3; %interictal

predClass = knnclassify(testFeats_nr,trainFeats,knnLabels,7);
predClass_rand = knnclassify(testFeats_rand,trainFeats_rand,knnLabels_rand,7);

accuracy = (sum(predClass == testLabels)/length(testLabels))*100
accuracyRand = (sum(predClass_rand == testLabels)/length(testLabels))*100


Sen = (sum(predClass == 1 | predClass == 2 & testLabels == 1 | testLabels == 2)/sum(testLabels == 1 | testLabels == 2))*100

Spec = (sum(predClass == 1 | predClass == 2 & testLabels == 3)/sum(testLabels == 3))*100
