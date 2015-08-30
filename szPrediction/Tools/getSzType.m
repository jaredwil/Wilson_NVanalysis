%Jared Wilson
%8/28/2015
% This is now a script that finds the sz annotations and return a matrix to
% be reference with the class of sz. Below is the legend for the type.
%    
%   0 - UCS (other electrical abnormality)
%   1 - CCS (clinically confirmed)
%   2 - CES (electricaly confirmed Sz)????
%   3 - NCS (not reported clinical Sz)
%
%%
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
szHorizon = 1; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%%
%begin function
%loop through all pts
for i = 12;  %%%%%TEMPORARY for debug

session = IEEGSession(pt{i},usernm,pswdBin) ;
fs = session.data.sampleRate;               %Find sampling Rate

timeUSec = cell(size(session.data.annLayer,2),1);
ptLabel = [];
startT = [];
%check to see if current pt. has annotations
if (size(session.data.annLayer,1) ~= 0)
    for layer = 1:size(session.data.annLayer,2)
        layerName = session.data.annLayer(layer).name;
        if(isempty(strfind(lower(layerName),'seizure')))
            continue;
        end

        %need to add to do this for only sz annot layers
        [~, timeUSec{layer}, ~] = getAllAnnots(session.data,layerName);
        
        if(isempty(strfind(lower(layerName),'ucs')) ~= 1)
            
            det = zeros(size(timeUSec{layer},1),1); %   0 - UCS (other electrical abnormality)
            ptLabel = [ptLabel; det];
            startT = [startT; (timeUSec{layer}(:,1)*1e-6)];
            
        elseif(isempty(strfind(lower(layerName),'ccs')) ~= 1)
            
            det = ones(size(timeUSec{layer},1),1);%   1 - CCS (clinically confirmed)
            ptLabel = [ptLabel; det];
            startT = [startT; (timeUSec{layer}(:,1)*1e-6)];    
            
        elseif(isempty(strfind(lower(layerName),'ces')) ~= 1)
            
            det = ones(size(timeUSec{layer},1),1)*2;%   2 - CES (electricaly confirmed Sz)????
            ptLabel = [ptLabel; det];
            startT = [startT; (timeUSec{layer}(:,1)*1e-6)];
            
        elseif(isempty(strfind(lower(layerName),'ncs')) ~= 1)
           
            det = ones(size(timeUSec{layer},1),1)*3;%   3 - NCS (not reported clinical Sz)
            ptLabel = [ptLabel; det];
            startT = [startT; (timeUSec{layer}(:,1)*1e-6)];

        end
    end

else
    %no annotations
    disp(['WARNING: NO Annotations for this Pt. ' ... 
        'startT/endT returned as an empty matrix'])
    startT = [];
    endT = [];
end

tot = [startT ptLabel];

%sort the start and end times from smalles to largest
[Y,Idx] = sort(tot(:,1));
sortTot = tot(Idx,:);



numTr = floor(length(startT)*trPct);  %number of training sz
%split the sz into train/test sets
trainLab = sortTot(1:numTr,:);
testLab = sortTot(numTr+1:end,:);


label = [pt{i} '_szPred_30secFeats.mat'];

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
saveLabel = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\30secFeats\' pt{i} '_szTLabels.mat'];
save(saveLabel,'szTLabels','-v7.3');




end