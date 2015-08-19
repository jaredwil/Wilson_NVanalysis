%Jared D. Wilson
%7/27/2015
%Regression pipeline Test -- Features Pre-Extracted to speed up testing
% ALL Features & Labels are loaded from .mat file contained on either
% dropbox or external hardrive on Laplace
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
addpath(genpath('H:\HumanNV\szPred_feats'))  %this is where .mat file are contained

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

    if(numWork < 2)  %processing is done on a laptop so don't do it in parallel
        parpool('local',1)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;
    elseif(numWork > 8) %limit the number of workers to 6
        parpool('local',8)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;
    else  %set up a parallel pool with max number of workers available between 2 and 6
        parpool(myCluster)
        p = gcp('nocreate'); % If no pool, do not create new one.
        poolsize = p.NumWorkers;

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%begin function
%loop through all pts0
for i = 12 %numel(pt)  %%%%%TEMPORARY for debug

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


%remove the interictal data so model is only trained on szhorizon for
%reression
% trainFeats(trainLabels > szHorizon*60*60,:) = [];
% trainLabels(trainLabels > szHorizon*60*60,:) = [];


%normalize features
avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);
trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);

%%
% Train Lasso Regression
disp(['Training Lasso Model on Patient: ' pt{i}])
opts = statset('UseParallel',true);  %Do dis shit in parallel!!!  use number of availble workers for CV
[fLasso, tINFO] = lasso(trainFeats, trainLabels,'CV',5,'Options',opts);
%BOOM!!!
disp(['DONE Training Lasso Model on Patient: ' pt{i}])

dzT = tINFO.DF; 
tInt = tINFO.Intercept;
% xAlpha = xINFO.Alpha;


%%
%%%%%%%%%%%%%%%%%%%%%%%  TEST   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%remove the interictal data so model is only trained on szhorizon for
%reression... train on all test on only szHorizon...
testFeats(testLabels > szHorizon*60*60,:) = [];
testLabels(testLabels > szHorizon*60*60,:) = [];
tInt2 = repmat(tInt, size(testFeats,1),1);

%normalize test feats
testFeats = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);

utLasso = testFeats*fLasso + tInt2;

tTest = repmat(testLabels,1,100);

tLasso_coor = corr(tTest,utLasso);
tLasso_coor = tLasso_coor(1,:);


minIdx = tINFO.IndexMinMSE;
seIdx  = tINFO.Index1SE;


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
name = [pt{i} '_lassoRes'];
% saveas(h,name,'jpg')


lassoPlot(fLasso,tINFO,'PlotType','CV');
% Use a log scale for MSE to see small MSE values better
% set(gca,'YScale','log');

figure(8)
% plot(utLasso(:,tLasso_coor == max(tLasso_coor)));
% plot(utLasso(:,minIdx));
plot(utLasso(:,seIdx));
hold on;
plot(testLabels);

% 
% bestLasso = fLasso(:,tLasso_coor == max(tLasso_coor));
% bestInt = tInt(tLasso_coor == max(tLasso_coor)); 
% numFeats = dzT(tLasso_coor == max(tLasso_coor)); 

% lassoRes = struct('coef',bestLasso,'int',bestInt,'numFeats',numFeats);
% saveLabel = [pt{i} '_bestLasso.mat'];
% save(saveLabel,'lassoRes','-v7.3');


end