%Time-2-Seizure Prediction Demov1
%Jared Wilson
%10/29/2015
%
% This demo script allows the user to select the pt. they wish to analyse
% from the NV dataset and uses existing LASSO models to produce results.
% The results displayed shows correlation between true labels and prediction for all
% seizures and and a plot comparing comparing the true labels to the
% predictions
%
% Feature windows contained in .mat files have all sz horizons.
%
% SzType labels contain information on the type of each sz:    
%   0 - UCS (other electrical abnormality)
%   1 - CCS (clinically confirmed)
%   2 - CES 
%   3 - NCS (not reported clinical Sz)
%
% ALL FILES Found in folder: H:\jaredwil\DEMO_files
% Models - T2Sz_models
% Data   - Win_5min

%PLEASE READ THROUGH ALL -- Some paramaters can be changed such as pt.
%number and are noted in code. Plots are provided to give user an
%understanding of what the model is doing. 

%%
% INIT
clear all; close all; clc;
warning('off')
addpath(genpath('Wilson_NVanalysis'))
addpath(genpath('H:\jaredwil\DEMO_files'))  %this is where .mat file are contained on local comp
set(0,'DefaultTextInterpreter','none');

szHorizon = 2; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THESE PARAMETERS CAN BE SET BY THE USER
i = 12;  %<- i CHOOSES WHICH PT. WILL BE EVALUATED
corrTh = .55;  %<- The seizures plotted will all have a correlation value greater than this threshold.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%START CODE
%text label for feat and type mat files
label = [pt{i} '_szPred_5minFeats.mat'];
tylabel = [pt{i} '_szTypeLabels5.mat'];

try   
    load(label);
    load(tylabel);
catch
    error('This pt. does not have a save .mat to load from... try a different pt.')
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% get sz time and type labels%%%%%%%%%%%%%%%%%%%%
trainSzT = szType.labels.train;
testSzT = szType.labels.test;

szTimes = [szType.szT.train;szType.szT.test];

szTimes(:,1) = szTimes(:,1)/60/60/24;

avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);

allFeat = [trainFeats; testFeats];
allFeat = bsxfun(@rdivide, bsxfun(@minus,allFeat,avgFeats), stdFeats);

allLabs  = [trainLabels;testLabels];
allSzLab  = [trainSzT;testSzT];


%%%%%%%%%%%%%% REMOVE SZ THAT ARE NOT CCS (CLINICALLY CONFIRMERD) %%%%%%%%%
%ONLY LOOK AT CCS
allFeat(allSzLab(:,2)~=1,:) = [];
allLabs(allSzLab(:,2)~=1) = [];
allSzLab(allSzLab(:,2)~=1,:) = [];

%load feature weights
featLab = ['H:\jaredwil\DEMO_files\T2Sz_models\' pt{i} '_bestLasso5.mat'];
load(featLab);  %feat label contains three models (MIN SSE, MAX CORR, AND 
%totmin (which is the lambda value 1 standard error from max corr, this 
%proved to be the best result I think because of regularization and avoided 
%overfitting)

% w - weights
% Int - offset (commonly called w_0)
w = lassoRes.coef.totmin;
tInt = lassoRes.int(3);
tInt2 = repmat(tInt, size(allFeat,1),1);

%prediction for all windows 
predLabOrg = allFeat*w + tInt2;


%The code below is used to find the correlation for each seizure
%individually. It is a little bit hard to follow and a bit sloppy. Unless
%interested skip.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 

szCorr = zeros(numel(numSz),1);
randCorr = zeros(numel(numSz),1000);

timeDay   = zeros(numel(numSz),1);
timeSec   = zeros(numel(numSz),1);
trainSz   = zeros(numel(numSz),1);
winDur_sz = zeros(numel(numSz),1);
for sz = 1:numel(numSz);

    winDur_sz(sz) = length(allLabs(uIdx == sz));
    szCorr(sz) = corr(allLabs(uIdx == sz),predLabOrg(uIdx == sz));

    timeDay(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
    timeSec(sz) = mode(allSzLab(uIdx == sz,1));

    trainSz(sz) = sum((timeSec(sz) == szType.szT.train(:,1)));
end

%eliminate sz that do not have many observations because corr results are
%not accurate 
idxElim = (winDur_sz < 6);
szCorr(idxElim)        = [];
timeDay(idxElim)       = [];
timeSec(idxElim)       = [];
trainSz(idxElim)       = [];

%PLOT SZ CORR FOR ALL SZ OVER ENTIRE DATASET
figure(1)
plot(timeDay(trainSz == 1),szCorr(trainSz == 1)*100,'b.','MarkerSize',20)
hold on;
plot(timeDay(trainSz == 0),szCorr(trainSz == 0)*100,'r.','MarkerSize',20)
vline(timeDay(trainSz == 1),'b:')
vline(timeDay(trainSz == 0),'r:')
grid on;
set(gca,'FontSize',12);
set(gca,'LineWidth',1);
set(gcf,'position',get(0,'screensize'));
vline(szType.szT.test(1)/60/60/24,'m-')
ylabel('Correlation (%)')
xlabel('Time (days)')
legend('Train','Test','location','best')
title('Correlation between True Time to Seizure Labels and Prediction for Each Seizure')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section creates a plot of predicted labels for test seizures based
% on the correlation threshold set at the top of the script
numSz = length(szCorr);

corrLab = zeros(size(allLabs));
isTest = zeros(size(allLabs));
retrainSz = zeros(numel(numSz),1);

trainSztimes = szType.szT.train(:,1);
trainSztimes = trainSztimes(trainSz == 1);

testSztimes = szType.szT.test(:,1);
% trainSztimes = trainSztimes(1:round(length(trainSztimes) * (0.4/0.7)));  %This reduces the training set from 70% to 30%

for sz = 1:numSz;
    corrLab(uIdx == sz) = szCorr(sz);                                   %repeat correlation for a label
    isTest(uIdx == sz) = (sum((timeSec(sz) == testSztimes)));  %sz was orignally trained on
%         retrainSz(sz) = sum((timeSec(sz) == trainSztimes)) & (szCorr(sz) > corrTh); %the retrain sz index is based on sz corr and time 
end

idxplotTest = (corrLab > corrTh) & isTest; 

plotFeats  = allFeat(idxplotTest,:);
plotTestLabels = allLabs(idxplotTest);
plotPred = predLabOrg(idxplotTest);

%visual assesment
time5 = linspace(5,length(plotPred)*5,length(plotPred))/60;
figure(2)
set(gca,'FontSize',15);
set(gca,'LineWidth',2);
set(gcf,'Position',get(0,'Screensize')); 
%make background white
set(gcf,'Color','w');
plot(time5, plotTestLabels/60/60);
hold on;
plot(time5, smooth(plotPred/60/60,5));
title(['Time to Seizure Visualization of Predictions with Corr. > ' num2str(corrTh*100) '%: pt. ' pt{i}])
xlabel('Time (Hours)')
ylabel('Time to Seizure (Hours)')
legend('Test','Prediction','location','best')





