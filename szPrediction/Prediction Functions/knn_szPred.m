function [ labels ] = knn_szPred(test_eeg)
%[ labels ] = knn_szPred(test_eeg, trainData, trainLabels)
%   This is a very simple prediction function but can be used as a model
%   for future predction functions. Inputs to this function should just be
%   the test_eeg and the model should be loaded in. 

addpath(genpath('C:\Users\Jared\Dropbox\Thesis stuff\szPrediction_models'))

load('trainFeats.mat')
load('trainLabels.mat')


labels = knnclassify(testFeats,trainFeats,trainLabels,20);

end

