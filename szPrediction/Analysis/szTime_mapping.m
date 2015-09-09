
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WARNING OLD DONT USE

%Jared Wilson
%8/25/15

%This script find the sz times for each pt. and maps all features contained
%in .mat files to corrosponding sz time to find sz cluster.

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
szHorizon = 2; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};


for i = 1:numel(pt)  %%%%%TEMPORARY for debug

%get sz start times
session = IEEGSession(pt{i},usernm,pswdBin) ;
fs = session.data.sampleRate;               %Find sampling Rate
 
%get ALL seziure start and end times from anotations
[startT, endT] = getSzTime(session);

if(isempty(startT))
   disp('NO ANNOTATION SKIP');
   continue;
end


numTr = floor(length(startT)*trPct);  %number of training sz
%split the sz into train/test sets
trainST = startT(1:numTr);
trainET = endT(1:numTr);
testST = startT(numTr+1:end);
testET = endT(numTr+1:end);    
    
label = [pt{i} '_szPred_5minFeats.mat'];

try   
    load(label);
catch
    disp('This pt. does not have a save .mat to load from')
    continue;

end

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

%combine labels
% labs = [trainLabels; testLabels];
trainszStIdx = [];
testszStIdx  = [];

[~,trainszStIdx] = findpeaks(trainLabels);
[~,testszStIdx] = findpeaks(testLabels);

%the first index is alwasy a start time of sz
trainszStIdx = [1; trainszStIdx];
testszStIdx = [1; testszStIdx];

szTrainLabels = zeros(size(trainLabels));
szTestLabels  = zeros(size(testLabels));
%create train labels
for sz = 1:length(trainszStIdx)
    if(sz == length(trainszStIdx))
        szTrainLabels(trainszStIdx(sz):end) = trainST(sz);
    else
        szTrainLabels(trainszStIdx(sz):trainszStIdx(sz+1)) = trainST(sz);
    end

end
%create test labels
for sz = 1:length(testszStIdx)
    if(sz == length(testszStIdx))
        szTestLabels(testszStIdx(sz):end) = testST(sz);
    else
        szTestLabels(testszStIdx(sz):testszStIdx(sz+1)) = testST(sz);
    end
end

szT = struct('train',trainST,'test',testST);
labels = struct('train',szTrainLabels,'test',szTestLabels);

szTLabels = struct('szT',szT,'labels',labels);
saveLabel = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\5minFeats\' pt{i} '_szTypeLabels5.mat'];
save(saveLabel,'szTLabels','-v7.3');

end

