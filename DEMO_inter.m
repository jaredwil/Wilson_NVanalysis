% This function allows for the model to be tested on interictal data
% Jared Wilson
% 10/29/2015
%
%
% WARNING: THIS FUNCTION DOES USE THE PORTAL SO IT IS GOING TO TAKE A VERY
% LONG TIME IF YOU WANT TO TEST ON ALOT OF DATA FOR EACH SEIZURE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INIT
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
addpath(genpath('H:\jaredwil\DEMO_files'))  %this is where .mat file are contained on local comp
set(0,'DefaultTextInterpreter','none');

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THESE ARE THE PARAMETERS TO CHANGE  (models trained on 2 hr hoizon)
extendedHorizon = 4; %This is the number of hours of data prior to each good seizure you want to test on.
i = 12;               %Pt. you want to test
corrTh = .55;         %The seizures that you will get data preceding will all have correlations above this threshold

% This is where you IEEG portal INFO goes
usernm = 'jaredwil';
pswdBin = 'jar_ieeglogin.bin';

szPlotDur = (extendedHorizon-1)*12-1; %number of prediction windows in plot 
%(this is used because some outtages could exist that could shrink the
%number of available windows) you might have to play with this value if you
%want to plot as much of the extended Horizon as possible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% begin function


    close all;

%   label = [pt{i} '_szPred_30secFeats.mat'];
    label = [pt{i} '_szPred_5minFeats.mat'];

    tylabel = [pt{i} '_szTypeLabels5.mat'];

    try   
        load(label);
        load(tylabel);
    catch
        disp('This pt. does not have a save .mat to load from')
        continue;

    end
    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate
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
    
    
%%%%%%%%%%%%%% REMOVE SZ THAT ARE NOT CCS (CLINICALLY CONFIRMERD) %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ONLY LOOK AT CCS
    allFeat(allSzLab(:,2)~=1,:) = [];
    allLabs(allSzLab(:,2)~=1) = [];
    allSzLab(allSzLab(:,2)~=1,:) = [];

    %load feature weights
    % load('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\resultsV1\NVC1001_23_005_bestLasso_CCS.mat')

%     featLab = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Res LassoLib\' pt{i} '_bestLasso5.mat'];
    featLab = [pt{i} '_bestLasso5.mat'];

    load(featLab);


    f = lassoRes.coef.totmin;
    tInt = lassoRes.int(3);
    tInt2 = repmat(tInt, size(allFeat,1),1);


    predLabOrg = allFeat*f + tInt2;

      %average feature across sz do pca then cluster. 
    [numSz,~, uIdx] = unique(allSzLab(:,1),'stable'); 
    
    szCorr = zeros(numel(numSz),1);
    randCorr = zeros(numel(numSz),1000);

    timeDay   = zeros(numel(numSz),1);
    timeSec   = zeros(numel(numSz),1);
    trainSz   = zeros(numel(numSz),1);
    winDur_sz = zeros(numel(numSz),1);
    for sz = 1:numel(numSz);
        
        winDur_sz(sz) = length(allLabs(uIdx == sz));
%       szAvg_feat(sz,:) = [mean(allFeat(uIdx == sz,:),1) szTimes(sz,1)];
        szCorr(sz) = corr(allLabs(uIdx == sz),predLabOrg(uIdx == sz));
        %random correlation

        
        timeDay(sz) = mode(allSzLab(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(allSzLab(uIdx == sz,1));
        
        trainSz(sz) = sum((timeSec(sz) == szType.szT.train(:,1)));
    end
    
    %eliminate sz that do not have many observationsbecause
    %it messes up the null hypothesis
    idxElim = (winDur_sz < 6);
    szCorrHist = szCorr;
    szCorrHist(idxElim)        = [];
    timeDay(idxElim)       = [];
%     timeSec(idxElim)       = [];
    trainSz(idxElim)       = [];
    
    
% Original corr by sz has been found now threshold and keep sz's with high
% correlation and plot results

    corrLab = zeros(size(allLabs));
    isTest = zeros(size(allLabs));
    retrainSz = zeros(numel(numSz),1);

    trainSztimes = szType.szT.train(:,1);
    trainSztimes = trainSztimes(trainSz == 1);
    
    testSztimes = szType.szT.test(:,1);
    goodSz_test = zeros(numel(numSz),1);
    for sz = 1:numel(numSz);
        corrLab(uIdx == sz) = szCorr(sz);                                   %repeat correlation for a label
        isTest(uIdx == sz) = (sum((timeSec(sz) == testSztimes)));  %sz was orignally trained on
        goodSz_test(sz) = timeSec(sz);
    end

    %these are the feature idx to train a new model on and it is important to
    %note this because they both are inside the new "training set" and are also
    %withing the correlation threshold for retraining
    idxplotTest = (corrLab > corrTh) & isTest; 

    plotFeats      = allFeat(idxplotTest,:);
    plotTestLabels = allLabs(idxplotTest);
    plotPred       = predLabOrg(idxplotTest);
    
    intTest_szLabels = allSzLab(idxplotTest,:);
    
    [szStartT] = unique(intTest_szLabels(:,1));
    szEndT     =  szStartT + 30;  %assume sz is about 30 seconds long this
                                    %is an assumption made cause figuring 
                                    %out exact time would take a long time
    winLen = 5*60;
    winDisp = 5*60;
    
    if(isempty(szStartT))
        error('NO GOOD SZ TO TEST ON');  %mosey along in the for loop if there are no good sz
    end
    [extFeats, ~] = getSzFeats(session, szStartT, szEndT, extendedHorizon, winLen, winDisp);
    %save the extFeats
    matName = ['H:\jaredwil\Lasso Results\szHorizon_test\' pt{i} '_szFeats.mat'];
    save(matName,'extFeats')
    
    
    testInt_lables = extFeats(:,1);
    
    testInt_feats =  extFeats(:,2:end);
    testInt_feats = bsxfun(@rdivide, bsxfun(@minus,testInt_feats,avgFeats), stdFeats); %NORMALIZE!!!
    tInt2 = repmat(tInt, size(testInt_feats,1),1);

    %make prediction and plot
    intPred = testInt_feats*f + tInt2;
    
    
    %visual assesment
    time5 = linspace(5,length(intPred)*5,length(intPred))/60;
    
    figure(1)
    set(gca,'FontSize',15);
    set(gca,'LineWidth',2);
    set(gcf,'Position',get(0,'Screensize')); 
    %make background white
    set(gcf,'Color','w');
    plot(time5, testInt_lables/60/60);
    hold on;
    plot(time5, smooth(intPred/60/60,5));
%     title(['Time to Seizure Visualization of Predictions with Better Than Random Corr.: pt. ' pt{i}])
    title(['Testing Model Beyond 2 Hour Trained Seizure Horizon: pt. ' pt{i}])

    xlabel('Time (Hours)')
    ylabel('Time to Seizure (Hours)')
    legend('Test','Prediction','location','best')
    plotName = ['H:\jaredwil\Lasso Results\szHorizon_test\' pt{i} '_plotResults'];
    saveas(gcf,plotName,'jpg')
    savefig(gcf,plotName)
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALL  THIS STUFF IS JUST FOR MAKING PLOTS FOR PREZENTATIONS!
    %plot each individual sz
    [~,testSzIdx] = findpeaks(testInt_lables);

    %the first index is alwasy a start time of sz
%     testSzIdx = [1; testSzIdx];
    szEndIdx = [testSzIdx-1];

%create train labels
szLabels = cell(length(szEndIdx+1),1);
predLabels = cell(length(szEndIdx+1),1);
for sz = 1:(length(szEndIdx)+1)
    if(sz == (length(szEndIdx)+1))
        szLabels{sz} = testInt_lables(end-szPlotDur:end);
        predLabels{sz} = intPred(end-szPlotDur:end);
    else
        szLabels{sz} = testInt_lables(szEndIdx(sz)-szPlotDur:szEndIdx(sz));
        predLabels{sz} = intPred(szEndIdx(sz)-szPlotDur:szEndIdx(sz));
    end

end
    
% THIS IS WHERE PLOTTING HAPPENS MAY WANT TO EDIT THIS STUFF DEPENDING ON
% WHAT KINDA PLOT YOU ARE LOOKING TO PRODUCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subplot
for sz = 1:(length(szEndIdx)+1)
    figure(10)
    
    subplot(length(szEndIdx)+1,1,sz)
    
    time5 = 0-(szLabels{sz}/60/60);
     set(gca,'FontSize',15);
    set(gca,'LineWidth',2);
    set(gcf,'Position',get(0,'Screensize')); 
    %make background white
    set(gcf,'Color','w');
    plot(time5, szLabels{sz}/60/60,'linewidth',2);
    hold on;
    plot(time5, smooth(predLabels{sz}/60/60,5),'linewidth',2);
    suptitle(['Testing Model Beyond 2 Hour Trained Seizure Horizon: pt. ' pt{i}])
    xlabel('Time (Hours)')
    ylabel('Time to Seizure (Hours)')
%     legend('Test','Prediction','location','best')
%     axis([-30 0 0 2])
    
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Individual Plots for each Sz

   for sz = 1:(length(szEndIdx)+1)
        figure(100+sz)

        time5 = 0-(szLabels{sz}/60/60);
         set(gca,'FontSize',15);
        set(gca,'LineWidth',2);
        set(gcf,'Position',get(0,'Screensize')); 
        %make background white
        set(gcf,'Color','w');
        plot(time5, szLabels{sz}/60/60);
        hold on;
        plot(time5, smooth(predLabels{sz}/60/60,5));
        %     title(['Time to Seizure Visualization of Predictions with Better Than Random Corr.: pt. ' pt{i}])
        title(['Testing Model Beyond 2 Hour Trained Seizure Horizon: pt. ' pt{i} ' Sz.' num2str(sz)])
        xlabel('Time (Hours)')
        ylabel('Time to Seizure (Hours)')
        legend('Test','Prediction','location','best')
        axis([min(time5) max(time5) 0 2])
    
    end 
    
    
    
