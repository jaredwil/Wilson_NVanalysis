%Jared D. Wilson
%7/27/2015
%Regression pipeline Test
%   This script is designed to creat a lasso regression model that will
%   serve two purposes. 
%   1.) Test regression for a seizure prediction problem
%   2.) Feature Selection

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
% trPct = 0.7;
trPct = 0.1;

winLen = 20;
winDisp = 20;
szHorizon = 1; %hours

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
%loop through all pts0
i = 12;  %%%%%TEMPORARY for debug

%Get train features and training labels (lables -> minutes to sz)
%features extracted defined in szPred_winFeatExt <- change this file
[trainFeats, trainLabels ] = szPred_train(pt{i}, usernm, pswdBin, trPct, winLen, winDisp, szHorizon);


%remove the interictal data so model is only trained on szhorizon for
%reression
trainFeats(trainLabels > szHorizon*60*60,:) = [];
trainLabels(trainLabels > szHorizon*60*60,:) = [];


%normalize features
avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);
trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);


%%
% Train Lasso Regression
[fLasso, tINFO] = lasso(trainFeats, trainLabels);

dzT = tINFO.DF; 
tInt = tINFO.Intercept;
% xAlpha = xINFO.Alpha;
tInt2 = repmat(tInt, size(trainFeats,1),1);


%%
%%%%%%%%%%%%%%%%%%%%%%%  TEST   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get the training features and labels
[testFeats, testLabels] = szPred_test(pt{i}, usernm, pswdBin, trPct, winLen, winDisp, szHorizon);

%remove the interictal data so model is only trained on szhorizon for
%reression
testFeats(testLabels > szHorizon*60*60,:) = [];
testLabels(testLabels > szHorizon*60*60,:) = [];

%normalize test feats
testFeats = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);

utLasso = testFeats*fLasso; + tInt2;

tTest = repmat(testLabels,1,100);

tLasso_coor = corr(tTest,utLasso);
tLasso_coor = tLasso_coor(1,:);

figure(6)
title('Number of Non-Zero Features vs Resulting Correlation')
plot(dzT,tLasso_coor*100,'o');
xlabel('Number of Non-Zero Feature Weights in Lasso Model')
ylabel('Resulting Test Correlation (%)')


bestLasso = fLasso(:,tLasso_coor == max(tLasso_coor));
bestInt = tInt(tLasso_coor == max(tLasso_coor)); 
