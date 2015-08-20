%Feature selection of classification scheme utilizing sequential feature
%selection
%Jared D. Wilson
%8/17/2015

%%
% Start
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
% addpath(genpath('H:\HumanNV\szPred_feats'))
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Lasso'))

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    elseif(numWork > 8) %limit the number of workers to 6
        parpool('local',8)
    else  %set up a parallel pool with max number of workers available between 2 and 6
        parpool(myCluster)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = cell(numel(pt),1);
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

    %normalize features
    avgFeats = mean(trainFeats,1);
    stdFeats = std(trainFeats,[],1);
    trainFeats = bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);

    %SVM Labels
    svmLabels = trainLabels;
    svmLabels(svmLabels >= 0 & svmLabels < szHorizon*60*60) = 1; %preictal
    svmLabels(svmLabels > szHorizon*60*60) = 0; %interictal

    % Train SVM model example   
%     model{i} = svmtrain(svmLabels, trainFeats, '-t 0 -b 1 -c 1'); 

    %%Create stuff to run an svm model
    %normalize test feats
    testFeats = bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);
    %SVM
    svmTestLabels = testLabels;
    svmTestLabels(svmTestLabels >= 0 & svmTestLabels < szHorizon*60*60) = 1; %preictal
    svmTestLabels(svmTestLabels > szHorizon*60*60) = 0; %interictal

    %now run sequentialfs 
    c = cvpartition([trainLabels;testLabels],'k',10);
    opts = statset('display','iter','UseParallel',true);
    fun = @(featTr,labTr,featTest,labTest)...   %train a new model based on subset of features selected
        (sum(~strcmp(labTest,svmpredict(labTest,featTest,svmtrain(labTr,featTr,'-t 0 -b 1 -c 1'),'-b 1'))));
        
    [fs,history] = sequentialfs(fun,[trainFeats; testFeats],[trainLabels;testLabels],'cv',c,'options',opts);

end