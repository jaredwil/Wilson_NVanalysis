%make better plots to show interictal results


clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\30secFeats'))  %this is where .mat file are contained on local comp
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\5minFeats'))  %this is where .mat file are contained on local comp
% addpath(genpath('H:\jaredwil\szPred_feats\Win_5min'))  %this is where .mat file are contained on local comp

set(0,'DefaultTextInterpreter','none');
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


for i = 1:numel(pt);  %%%%%TEMPORARY for debug
    close all;
             weightLabel = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Res LassoLib\' pt{i} '_bestLasso5.mat'];
        %     weightLabel = ['H:\jaredwil\Lasso Results\Res 8-29\' pt{i} '_bestLasso5.mat'];
             featLabel = ['C:\Users\Jared\Dropbox\NVanalysis_data\extendedSzHorizon_Test\' pt{i} '_szFeats'];
    try
        load(weightLabel);
        load(featLabel);
    catch
        disp('no data available here')
        continue;
    end
    
    if(isempty(extFeats))
       continue; 
    end
    
    
     label = [pt{i} '_szPred_5minFeats.mat'];

    tylabel = [pt{i} '_szTypeLabels5.mat'];

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
    
    avgFeats = mean(trainFeats,1);
    stdFeats = std(trainFeats,[],1);

    labels = extFeats(:,1);
    feats = extFeats(:,2:end);
    feats = bsxfun(@rdivide, bsxfun(@minus,feats,avgFeats), stdFeats);

    f = lassoRes.coef.totmin;
    tInt = lassoRes.int(3);
    tInt2 = repmat(tInt, size(feats,1),1);

    pred = feats*f + tInt2;

    figure(100)
    plot(labels)
    hold on;
    plot(pred)

    [PKS,LOCS]= findpeaks(labels);

    stSz  = [1; (LOCS)];
    endSz = [LOCS-1; length(labels)];

    trueSz = cell(numel(stSz),1);
    predSz = cell(numel(stSz),1);
    for sz = 1:numel(stSz)
       tmpTrue = labels(stSz(sz):endSz(sz))/60/60;
       tmpPred = pred(stSz(sz):endSz(sz))/60/60;
       xLab = flipud(tmpTrue);

       figure(sz)
        set(gca,'LineWidth',2);
        set(gcf,'Position',get(0,'Screensize')); 
        %make background white
        set(gcf,'Color','w');
       plot(xLab,tmpTrue);
       hold on;
       grid on;
       plot(xLab,tmpPred);
       legend('True Label','Prediction','location','best')
       title(['Time to Seizure Prediction: pt. ' pt{i} ' sz.' num2str(sz)])

        xlabel('Time (Hours)')
        ylabel('Time to Seizure (Hours)')
        legend('Test','Prediction','location','best')
        set(gca,'FontSize',15);
        axis([0 6-(5/60) 0 6-(5/60)])

        plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\extendedSzHorizon_Test\' pt{i} '\predRes_sz' num2str(sz)];
        saveas(gcf,plotName,'jpg')
%         savefig(gcf,plotName)

    end
        

        
end