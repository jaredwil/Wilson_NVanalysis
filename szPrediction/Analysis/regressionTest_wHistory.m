%Jared D. Wilson
%7/27/2015
%Regression pipeline Test -- Features Pre-Extracted to speed up testing
% ALL Features & Labels are loaded from .mat file contained on either
% dropbox or external hardrive on Laplace
%   This script is designed to creat a lasso regression model that will
%   serve two purposes. 
%   1.) Test regression for a seizure prediction problem
%   2.) Feature Selection
%
%NOTE: This script trains only on defined preictal data (2 hours before sz
%start time), Interictal features are not contained in script
%NVC...5minFeats but are contained int NVC...30secFeats.

%%
% Start
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
addpath(genpath('H:\jaredwil\szPred_feats\Win_5min'))  %this is where .mat file are contained on Laplace
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Lasso'))  %this is where .mat file are contained on local comp

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


lassoTime = zeros(1,numel(pt));
%%
%begin function
%loop through all pts0
for i = 1:numel(pt)  %%%%%TEMPORARY for debug
close all;
clear ('h','h1','h2');
    
    
label = [pt{i} '_szPred_5minFeats.mat'];
tylabel = [pt{i} '_szTypeLabels5.mat'];

try   
    load(label);
    load(tylabel);

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

%Remove intericatal Data from training data
trainFeats(trainLabels > szHorizon*60*60,:) = [];
trainLabels(trainLabels > szHorizon*60*60,:) = [];

%normalize features
avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);
trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);

%%%%%%%%%%%%%%%%%%%%%%%%%%% get sz time and type labels%%%%%%%%%%%%%%%%%%%%

    trainSzT = szType.labels.train;
    testSzT = szType.labels.test;
    
    szTimes = [szType.szT.train;szType.szT.test];
    
    szTimes(:,1) = szTimes(:,1)/60/60/24;
    
    allFeat = [trainFeats; testFeats];
    allSzLab  = [trainSzT;testSzT];
    
    
%%%%%%%%%%%%%% REMOVE SOME SZ %%%%%%%%%%%%%%%%%%%%%%%%%%%
%ONLY LOOK AT CCS
trainFeats(trainSzT(:,2)~=1,:) = [];
trainLabels(trainSzT(:,2)~=1,:) = [];
testFeats(testSzT(:,2)~=1,:) = [];
testLabels(testSzT(:,2)~=1,:) = [];

% %Remove UCS
% trainFeats(trainSzT(:,2)==0,:) = [];
% trainLabels(trainSzT(:,2)==0,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%RESHAPE TRAIN FEATURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1:2; %number of history samples
samples = size(trainFeats,1);
features = size(trainFeats,2);
startSAMP = max(N)+1;
M = (samples-length(N));  %number of time bins
%create 'R' matrix for linear regression algorithm
r = zeros(M, features*length(N)+1);
for resIdx = 1:M
    temp = trainFeats(startSAMP + (resIdx-1) - N,:);   %temp is a temporary matrix    
    r(resIdx,:) = [1 temp(:)'];
end

trainLabels = trainLabels(startSAMP:end);

[~,trainIdxRemove] = findpeaks(trainLabels);
[~,testIdxRemove] = findpeaks(testLabels);

r(trainIdxRemove,:) = [];
trainLabels(trainIdxRemove) = [];

trainR = r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.75;
%%
% Train Lasso Regression UseParallel
tic;
disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
opts = statset('UseParallel',true);  %Do dis shit in parallel!!!  use number of availble workers for CV
% [fLasso, tINFO] = lasso(trainFeats, trainLabels,'CV',5,'Options',opts);
[fLasso, tINFO] = lasso(trainR, trainLabels,'Alpha',alpha,'CV',3,'NumLambda',25,'Options',opts);
%BOOM!!!
disp(['DONE Training Lasso Model on Patient: ' pt{i}])
lassoTime(i) = toc;

dzT = tINFO.DF; 
tInt = tINFO.Intercept;
% xAlpha = xINFO.Alpha;
tInt2 = repmat(tInt, size(trainR,1),1);

utLasso = trainR*fLasso + tInt2;

tTrain = repmat(trainLabels,1,100);

tLasso_coor = corr(tTrain,utLasso);
tLasso_coor = tLasso_coor(1,:);

minIdx = tINFO.IndexMinMSE;
seIdx  = tINFO.Index1SE;
corrIdx = tLasso_coor == max(tLasso_coor);

figure(6)
title('Number of Non-Zero Features vs Resulting Correlation')
h = plot(dzT,tLasso_coor*100,'r.');
xlabel('Number of Non-Zero Feature Weights in Lasso Model')
ylabel('Resulting Test Correlation (%)')
grid on;
hold on;
h = plot(dzT(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
h = vline(dzT(minIdx),'g:');
h = plot(dzT(seIdx),tLasso_coor(seIdx)*100,'bs','LineWidth',4);
h = vline(dzT(seIdx),'b:');


h2 = lassoPlot(fLasso,tINFO,'PlotType','CV');
% Use a log scale for MSE to see small MSE values better
% set(gca,'YScale','log');
plotName = ['H:\jaredwil\Lasso Results\matLasso_history\' pt{i} '_mseRes_CCS'];
saveas(h2,plotName,'jpg')


bestLasso_corr = fLasso(:,corrIdx);
bestLasso_1SE = fLasso(:,seIdx);
bestLasso_Min = fLasso(:,minIdx);

bestInt_corr = tInt(corrIdx); 
bestInt_1SE = tInt(seIdx); 
bestInt_Min = tInt(minIdx); 

numFeats_corr = dzT(corrIdx); 
numFeats_1SE = dzT(seIdx); 
numFeats_Min = dzT(minIdx);

%SAVE ALL 
lassoRes = struct('coef',struct('corr',bestLasso_corr,'min',bestLasso_Min, ...
    'se',bestLasso_1SE),'int',[bestInt_corr bestInt_Min bestInt_1SE], ...
    'numFeats', [numFeats_corr numFeats_Min numFeats_1SE]);
saveLabel = ['H:\jaredwil\Lasso Results\matLasso_history\' pt{i} '_bestLasso_CCS.mat'];
save(saveLabel,'lassoRes','-v7.3');


%%
%%%%%%%%%%%%%%%%%%%%%%%  TEST   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%normalize test feats
testFeats = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);
%%%%%%%%%%%%%%%%%%%%%%%RESHAPE TEST FEATURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
samples = size(testFeats,1);
features = size(testFeats,2);
startSAMP = max(N)+1;
M = (samples-length(N));  %number of time bins
%create 'R' matrix for linear regression algorithm
r = zeros(M, features*length(N)+1);
temp = 1;
for resIdx = 1:M
    temp = testFeats(startSAMP + (resIdx-1) - N,:);   %temp is a temporary matrix    
    r(resIdx,:) = [1 temp(:)'];
    temp = 1;
end

testLabels = testLabels(startSAMP:end);
testR = r;
tInt2 = repmat(tInt, size(testR,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

utLasso = testR*fLasso + tInt2;


tTest = repmat(testLabels,1,size(utLasso,2));

tLasso_coor = corr(tTest,utLasso);
tLasso_coor = tLasso_coor(1,:);

% tLasso_MSE = mean(bsxfun(@minus,utLasso,testLabels).^2,1);
% tLasso_MSE = mean(abs(bsxfun(@minus,utLasso,testLabels)),1);

% 
% figure(10)
% plot(lambda, tLasso_MSE,'r.-')

% minIdx = tINFO.IndexMinMSE;
% seIdx  = tINFO.Index1SE;

figure(7)
title('Number of Non-Zero Features vs Resulting Correlation')
h = plot(dzT,tLasso_coor*100,'b.-');
% legend('Train','Test')
xlabel('Number of Non-Zero Feature Weights in Lasso Model')
ylabel('Resulting Test Correlation (%)')
grid on;
hold on;
h = plot(dzT(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
h = vline(dzT(minIdx),'g:');
h = plot(dzT(seIdx),tLasso_coor(seIdx)*100,'bs','LineWidth',4);
h = vline(dzT(seIdx),'b:');
plotName = ['H:\jaredwil\Lasso Results\matLasso_history\' pt{i} '_corrRes_CCS'];
saveas(h,plotName,'jpg')

uSE = testR*fLasso(:,seIdx) + tInt2(seIdx);
uMin = testR*fLasso(:,minIdx) + tInt2(minIdx);
uCorr = testR*fLasso(:,corrIdx) + tInt2(corrIdx);

%TEST/PLOT RESULTS!
timeMin = linspace(0,(size(testLabels,1)-1)*30,size(testLabels,1))/60;

figure(9)
h3 = plot(timeMin,smooth(uSE/60,5));
hold on;
h3 = plot(timeMin,smooth(uMin/60,5));
h3 = plot(timeMin,smooth(uCorr/60,5));
h3 = plot(timeMin,testLabels/60);
hline(30,'r:');
ylabel('Time to Sz (min)')
xlabel('Time (min)')
xlim([0 max(timeMin)])
legend('SE','Min','Corr','TEST')
% ylim([min(testLabels/60)-std(testLabels/60) , max(testLabels/60)+std(testLabels/60)])
plotName = ['H:\jaredwil\Lasso Results\matLasso_history\' pt{i} '_lassoRes_CCS'];
saveas(h3,plotName,'fig')
saveas(h3,plotName,'jpg')


end