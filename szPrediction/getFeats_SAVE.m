%GET ALL DEM FEATURES AND SAVE EM!!!!!
%Jared Wilson
%8/12/15

%%
% Start
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

% winLen = 30;
% winDisp = 30;
% szHorizon = 2; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%%
%begin function
%loop through all pts.
for i = 11%:numel(pt);  %%%%%TEMPORARY for debug

%Get train features and training labels (lables -> minutes to sz)
[trainFeats, trainLabels ] = szPred_train(pt{i}, usernm, pswdBin, trPct, winLen, winDisp, szHorizon);

%This means there are no labels so skip this pt.
if(isempty(trainFeats))
    continue;
end

%normalize features
avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);

%%
%%%%%%%%%%%%%%%%%%%%%%%  TEST   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get the training features and labels
[testFeats, testLabels] = szPred_test(pt{i}, usernm, pswdBin, trPct, winLen, winDisp, szHorizon);

data = struct('train',[trainLabels trainFeats],'mean',avgFeats,'std',stdFeats,'test',[testLabels testFeats]);
saveLabel = [pt{i} '_szPred_30secFeats.mat'];
save(saveLabel,'data','-v7.3');

end
