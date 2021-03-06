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
addpath(genpath('H:\jaredwil\szPred_feats'))  %this is where .mat file are contained on Laplace
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp

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

%%%%%%%%%%%%%%%%%%%%%% start par pool%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Establish parpool
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)        %if not p will contain all info about current pool;
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
%initialize paralell pool if available
if(poolsize == 0)

    myCluster = parcluster('local');
    numWork = myCluster.NumWorkers;

    if(numWork <= 2)  %processing is done on a laptop so don't do it in parallel
        parpool('local',1)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;
    elseif(numWork > 10) %limit the number of workers to 6
        parpool('local',10)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;
    else  %set up a parallel pool with max number of workers available between 2 and 6
        parpool(myCluster)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lassoTime = zeros(1,numel(pt));
%%
%begin function
%loop through all pts0
for i = 12%:numel(pt)  %%%%%TEMPORARY for debug
close all;
clear ('h','h1','h2');
    
    
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

%Remove intericatal Data from training data
trainFeats(trainLabels > szHorizon*60*60,:) = [];
trainLabels(trainLabels > szHorizon*60*60,:) = [];

testFeats(testLabels > szHorizon*60*60,:) = [];
testLabels(testLabels > szHorizon*60*60,:) = [];
%normalize features
avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);
trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);

%%%%%%%%%%%%%%%%%%%%%%%RESHAPE TRAIN FEATURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1:4; %number of history samples
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

trainR = r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numFolds = 10;
cvIdx    = crossvalind('Kfold', size(trainR,1), numFolds);
numLam   = 100;
lambda   = linspace(1e4,1e6,numLam);

tLasso_coor = zeros(numFolds,numLam);
tLasso_MSE  = zeros(numFolds,numLam);
%cross validation to find the best lambda
for cvIter = 1:numFolds
    disp(['Cross-Validation Progress: ' num2str(cvIter) '/' num2str(numFolds)])
    trainFeatsCV = trainR(cvIdx ~= cvIter,:);
    trainLabelsCV = trainLabels(cvIdx ~= cvIter,:);
    
    evalFeats = trainR(cvIdx == cvIter,:);
    evalLabels = trainLabels(cvIdx == cvIter,:);
    %%
    % Train Lasso Regression UseParallel
    w        = cell(1,numLam);
    numIter  = cell(numLam,1);

    tic;
    disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
    parfor_progress(numLam);
    parfor lam = 1:length(lambda)

    %     disp(['Progress: '  num2str(lam) '/' num2str(length(lambda))]);
        [w{lam}, ~, numIter{lam}] = LassoBlockCoordinate(trainFeatsCV,trainLabelsCV,lambda(lam),'maxIter',50000);
    %     [w{lam}, wp, numIter{lam}] = LassoIteratedRidge(trainFeats,trainLabels,lambda(lam),'maxIter',50000);
        parfor_progress;

    end
    parfor_progress(0);

    %BOOM!!! Done.
    disp(['DONE Training Lasso Model on Patient: ' pt{i}])
    lassoTime(i) = toc;

    w = cell2mat(w);

    dzT = sum(w == 0,1);
    %compute intercepts
    tInt = repmat(mean(trainLabelsCV),1,numLam) - mean(trainFeatsCV*w,1);

    tInt2 = repmat(tInt, size(evalFeats,1),1);

    utLasso = evalFeats*w + tInt2;

    tTrain = repmat(evalLabels,1,numLam);

    tmp = corr(tTrain,utLasso);
    tLasso_coor(cvIter,:) = tmp(1,:);

    tLasso_MSE(cvIter,:) = mean(bsxfun(@minus,utLasso,evalLabels).^2,1);

end
stdMSE = std(tLasso_MSE);
tLasso_MSE = mean(tLasso_MSE,1);

tLasso_coor = mean(tLasso_coor,1);

figure(1)
semilogx(lambda, tLasso_MSE,'r.-')
xlabel('Lambda')
ylabel('MSE')


% seIdx  = tINFO.Index1SE;
corrIdx = tLasso_coor == max(tLasso_coor);
minIdx = tLasso_MSE == min(tLasso_MSE);

figure(6)
title('Number of Non-Zero Features vs Resulting Correlation')
h = semilogx(lambda,tLasso_coor*100,'r.-');
% xlabel('Number of Feature Weights Zeroed in Lasso Model')
xlabel('Lambda')
ylabel('Resulting Test Correlation (%)')
grid on;
hold on;
% h = plot(dzT(minIdx),tLasso_coor(minIdx)*100,'gs','LineWidth',4);
% h = vline(dzT(minIdx),'g:');
% h = plot(dzT(seIdx),tLasso_coor(seIdx)*100,'bs','LineWidth',4);
% h = vline(dzT(seIdx),'b:');
plotName = ['H:\jaredwil\Lasso Results\' pt{i} '_corrRes'];
saveas(h,plotName,'jpg')



% h2 = lassoPlot(fLasso,tINFO,'PlotType','CV');
% Use a log scale for MSE to see small MSE values better
% set(gca,'YScale','log');
% plotName = ['H:\jaredwil\Lasso Results\' pt{i} '_mseRes'];
% saveas(h2,plotName,'jpg')

%Retrain with best lambda
corrLam = lambda(corrIdx);
minLam = lambda(minIdx);
        
[bestLasso_corr, ~, ~] = LassoBlockCoordinate(trainR,trainLabels,lambda(corrIdx),'maxIter',50000);
[bestLasso_Min, ~, ~] = LassoBlockCoordinate(trainR,trainLabels,lambda(minIdx),'maxIter',50000);


% bestLasso_corr = w(:,corrIdx);
% % bestLasso_1SE = fLasso(:,seIdx);
% bestLasso_Min = w(:,minIdx);

numFeats_corr = sum(bestLasso_corr == 0);
% numFeats_1SE = dzT(seIdx); 
numFeats_Min = sum(bestLasso_Min == 0);


int_corr = mean(trainLabels) -  mean(trainR*bestLasso_corr);
int_Min  = mean(trainLabels) -  mean(trainR*bestLasso_Min);


% figure(55)
% plot(trainLabels);
% hold on;
% plot(utLasso(:,minIdx));


%SAVE ALL 
lassoRes = struct('coef',struct('corr',bestLasso_corr,'min',bestLasso_Min), ...
    'int',[int_corr int_Min], ...
    'numFeats', [numFeats_corr numFeats_Min]);
saveLabel = ['H:\jaredwil\Lasso Results\' pt{i} '_bestLasso.mat'];
save(saveLabel,'lassoRes','-v7.3');


%%
%%%%%%%%%%%%%%%%%%%%%%%  TEST   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tInt2 = repmat(tInt, size(testFeats,1),1);

%normalize test feats
testFeats = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);
%%%%%%%%%%%%%%%%%%%%%%%RESHAPE TRAIN FEATURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

testR = r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tInt2 = repmat(tInt, size(testR,1),1);

%TEST SET CORR TEST
utLasso = testR*w + tInt2;

tTest = repmat(testLabels,1,numLam);

tLasso_coor = corr(tTest,utLasso);
tLasso_coor = tLasso_coor(1,:);

tLasso_MSE = mean(bsxfun(@minus,utLasso,testLabels).^2,1);

figure(10)
semilogx(lambda, tLasso_MSE,'r.-')
xlabel('Lambda')
ylabel('MSE')
% minIdx = tINFO.IndexMinMSE;
% seIdx  = tINFO.Index1SE;

figure(6)
hold on;
title('Number of Non-Zero Features vs Resulting Correlation')
h = semilogx(lambda,tLasso_coor*100,'b.-');
legend('Train','Test')


% uSE = testFeats*fLasso(:,seIdx) + tInt2(seIdx);
uMin = testR*w(:,minIdx) + tInt2(minIdx);
uCorr = testR*w(:,corrIdx) + tInt2(:,corrIdx);
% uCorr = testFeats*w(:,5);

%TEST/PLOT RESULTS!
timeMin = linspace(0,(size(testLabels,1)-1)*30,size(testLabels,1))/60;

figure(9)
% h3 = plot(timeMin,uSE/60);
h3 = plot(timeMin,uMin/60);
hold on;
h3 = plot(timeMin,smooth(uCorr/60,10));
h3 = plot(timeMin,testLabels/60);
hline(30,'r:');
ylabel('Time to Sz (min)')
xlabel('Time (min)')
xlim([0 max(timeMin)])
ylim([min(testLabels/60)-std(testLabels/60) , max(testLabels/60)+std(testLabels/60)])
legend('Min MSE','Max Correlation','Test Labels')
plotName = ['H:\jaredwil\Lasso Results\' pt{i} '_lassoRes'];
saveas(h3,plotName,'jpg')

end