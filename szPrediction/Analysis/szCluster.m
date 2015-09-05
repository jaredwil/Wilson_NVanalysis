%Jared Wilson
%8/30/2015

%clust seizure using k-means

%%
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\30secFeats'))  %this is where .mat file are contained on local comp
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\5minFeats'))  %this is where .mat file are contained on local comp


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
for i = 7;  %%%%%TEMPORARY for debug
       
%     label = [pt{i} '_szPred_30secFeats.mat'];
    label = [pt{i} '_szPred_5minFeats.mat'];

    tylabel = [pt{i} '_szTypeLabels.mat'];

    try   
        load(label);
        load(tylabel);
    catch
        disp('This pt. does not have a save .mat to load from')
        continue;

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

    %Remove intericatal Data from training data & Testing Data
    trainFeats(trainLabels > szHorizon*60*60,:) = [];
    trainLabels(trainLabels > szHorizon*60*60,:) = [];
    testFeats(testLabels > szHorizon*60*60,:) = [];
    testLabels(testLabels > szHorizon*60*60,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%% get sz time and type labels%%%%%%%%%%%%%%%%%%%%

    trainSzT = szType.labels.train;
    testSzT = szType.labels.test;
    
    szTimes = [szType.szT.train;szType.szT.test];
    
    szTimes(:,1) = szTimes(:,1)/60/60/24;
    
    allFeat = [trainFeats; testFeats];
    allSzLab  = [trainSzT;testSzT];
    
    
%%%%%%%%%%%%%% REMOVE SOME SZ %%%%%%%%%%%%%%%%%%%%%%%%%%%
%ONLY LOOK AT CCS
allFeat(allSzLab(:,2)~=1,:) = [];
allSzLab(allSzLab(:,2)~=1,:) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %szTimeline
    figure(1)
    vline(szTimes(:,1))
    
    %average feature across sz do pca then cluster. 
    [numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 
    
%     szAvg_feat = zeros(numel(numSz),(size(trainFeats,2) + 1));  %add time as a feature
    szAvg_feat = zeros(numel(numSz),size(trainFeats,2));  %add time as a feature

    for sz = 1:numel(numSz);
       
%         szAvg_feat(sz,:) = [mean(allFeat(uIdx == sz,:),1) szTimes(sz,1)];
        szAvg_feat(sz,:) = mean(allFeat(uIdx == sz,:),1);

    end
    
    %do some PCA shiz
    %do pca anlysis
    [coef, score, latent] = pca(szAvg_feat);

    %use top two principal components
    pcOne = score(:,1);
    pcTwo = score(:,2);
    pcThr = score(:,3);
    pcFou = score(:,4);
    pcFve = score(:,5);
    pcSix = score(:,6);

    %plot two components to show seperations
    figure(8)
    scatter(pcOne,pcTwo);
    xlabel('Principal Component One')
    ylabel('Principal Component Two')
    title('Spike Waveforms Scatter Plot Represented by Top Two Principal Components')
    
    %Show explained variance as a function of each principal component
    figure(9)
    plot((latent./sum(latent(:)))*100, 'bo')
    xlim([0 63])
    xlabel('Principal Component')
    ylabel('Variance Explained (%)')
    title('Principal Component vs. Total Variance Explained')
    explVar = sum(latent(1:2))/sum(latent(:))

    
    pcOne = (pcOne - mean(pcOne))./std(pcOne);
    pcTwo = (pcTwo - mean(pcTwo))./std(pcTwo);
    pcThr = (pcThr - mean(pcThr))./std(pcThr);
    pcFou = (pcFou - mean(pcFou))./std(pcFou);
    pcFve = (pcFve - mean(pcFve))./std(pcFve);
    pcSix = (pcSix - mean(pcSix))./std(pcSix);

    pc = [];
    pc = [pcOne pcTwo pcThr];% pcFou pcFve pcSix];

    %Clustering of spikes using k means and top two pc's
    group = kmeans(pc,2,'Distance','correlation','Replicates',100);
    
    
        
    %szTimeline v2 with group
    figure(11)
    vline(szTimes(group == 1,1),'r:')
    vline(szTimes(group == 2,1),'b:')
%     vline(szTimes(group == 3,1),'g:')

    
    
%     Not sure this really shows anything    
    figure(10)
    gscatter(pcOne,pcTwo,group,'rgbk')
    legend('Spike Group 1','Spike Group 2','Location','best')
    title('Spike Waveforms Clustered by Top Two Principal Components')

     
    
    
    




end