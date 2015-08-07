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
winLen = 5;
winDisp = 2.5;
szHorizon = 2; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%%
%begin function
%loop through all pts
i = 12;  %%%%%TEMPORARY for debug

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

%%
%HERE IS WHERE YOU WOULD TRAIN YOUR MODEL
%simple K-NN classification classifying 2hr (2) vs 1hr (1) vs
%interictal (3)

% knnLabels = trainLabels;
% knnLabels(knnLabels >= 0 & knnLabels < 1*60*60) = 1; %between 2-1 hours
% knnLabels(knnLabels >= 1*60*60 & knnLabels < szHorizon*60*60) = 2; %between 2-1 hours
% knnLabels(knnLabels > szHorizon*60*60) = 3; %interictal


% Train SVM model

model = svmtrain(svmLabels, trainFeats, '-t 0 -b 1 -c 1'); 


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

[predLab, accuracy, prob_estimates] = svmpredict(svmTestLabels, testFeats, model, '-b 1');

%simple test K-nn
%CREATE LABELS  classification classifying 1hr (1) vs 2hr (2) vs interictal (3)
% testLabels(testLabels >= 0 & testLabels < 1*60*60) = 1; %between 2-1 hours
% testLabels(testLabels >= 1*60*60 & testLabels < szHorizon*60*60) = 2; %between 2-1 hours
% testLabels(testLabels > szHorizon*60*60) = 3; %interictal
% 
% predClass = knnclassify(testFeats,trainFeats,knnLabels,7);
% 
% accuracy = (sum(predClass == testLabels)/length(testLabels))*100
% 
% 
% Sen = (sum(predClass == 1 | predClass == 2 & testLabels == 1 | testLabels == 2)/sum(testLabels == 1 | testLabels == 2))*100
% 
% Spec = (sum(predClass == 1 | predClass == 2 & testLabels == 3)/sum(testLabels == 3))*100
%  
% 

