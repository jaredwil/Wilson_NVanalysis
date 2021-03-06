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

%%%%%%%%%%%%%%%%%%%%%% start par pool%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Establish parpool
% p = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(p)        %if not p will contain all info about current pool;
%     poolsize = 0;
% else
%     poolsize = p.NumWorkers;
% end
% %initialize paralell pool if available
% if(poolsize == 0)
% 
%     myCluster = parcluster('local');
%     numWork = myCluster.NumWorkers;
% 
%     if(numWork <= 2)  %processing is done on a laptop so don't do it in parallel
%         parpool('local',1)
%         p = gcp('nocreate'); % If no pool, do not create new one.
%         poolsize = p.NumWorkers;
%     elseif(numWork > 6) %limit the number of workers to 6
%         parpool('local',6)
%         p = gcp('nocreate'); % If no pool, do not create new one.
%         poolsize = p.NumWorkers;
%     else  %set up a parallel pool with max number of workers available between 2 and 6
%         parpool(myCluster)
%         p = gcp('nocreate'); % If no pool, do not create new one.
%         poolsize = p.NumWorkers;
% 
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lassoTime = zeros(1,numel(pt));
%%
%begin function
%loop through all pts0
for i = 12 %numel(pt)  %%%%%TEMPORARY for debug
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

%normalize features
avgFeats = mean(trainFeats,1);
stdFeats = std(trainFeats,[],1);
trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);

%%
% Train Lasso Regression UseParallel
tic;
disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
opts = statset('UseParallel',true);  %Do dis shit in parallel!!!  use number of availble workers for CV
% [fLasso, tINFO] = lasso(trainFeats, trainLabels,'CV',5,'Options',opts);
[fLasso, tINFO] = lasso(trainFeats, trainLabels,'CV',2,'Options',opts);
%BOOM!!!
disp(['DONE Training Lasso Model on Patient: ' pt{i}])
lassoTime(i) = toc;

dzT = tINFO.DF; 
tInt = tINFO.Intercept;
% xAlpha = xINFO.Alpha;
tInt2 = repmat(tInt, size(trainFeats,1),1);

utLasso = trainFeats*fLasso + tInt2;

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
plotName = ['H:\jaredwil\Lasso Results' pt{i} '_corrRes_regLasso'];
saveas(h,plotName,'jpg')

% h2 = lassoPlot(fLasso,tINFO,'PlotType','CV');
% % Use a log scale for MSE to see small MSE values better
% % set(gca,'YScale','log');
% plotName = ['H:\jaredwil\Lasso Results' pt{i} '_mseRes'];
% saveas(h2,plotName,'jpg')


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
saveLabel = ['H:\jaredwil\Lasso Results' pt{i} '_bestLassoRegLasso.mat'];
save(saveLabel,'lassoRes','-v7.3');


%%
%%%%%%%%%%%%%%%%%%%%%%%  TEST   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tInt2 = repmat(tInt, size(testFeats,1),1);

%normalize test feats
testFeats = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);

uSE = testFeats*fLasso(:,seIdx) + tInt2(seIdx);
uMin = testFeats*fLasso(:,minIdx) + tInt2(minIdx);
uCorr = testFeats*fLasso(:,corrIdx) + tInt2(corrIdx);

%TEST/PLOT RESULTS!
timeMin = linspace(0,(size(testLabels,1)-1)*30,size(testLabels,1))/60;

figure(9)
h3 = plot(timeMin,smooth(uSE/60,5));
hold on;
h3 = plot(timeMin,smooth(uMin/60,5));
% h3 = plot(timeMin,uCorr/60);
h3 = plot(timeMin,testLabels/60);
hline(30,'r:');
ylabel('Time to Sz (min)')
xlabel('Time (min)')
xlim([0 max(timeMin)])
legend('SE','Min','Corr','TEST')
% ylim([min(testLabels/60)-std(testLabels/60) , max(testLabels/60)+std(testLabels/60)])
plotName = ['H:\jaredwil\Lasso Results' pt{i} '_lassoRes_regLasso'];
saveas(h3,plotName,'jpg')




end