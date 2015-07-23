%Jared D. Wilson
% 7/17/2015
%This script analyzes the hump and spikes that are present in various
%features extracted from the NV dataset.

% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
%add path to all mat file with feature file in them.
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\allCh_2months_OneminWinFeats'))
set(0, 'DefaulttextInterpreter', 'none') 

pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%define all time in seconds
day = 86400; %sec
hour = 3600; %sec
min = 60; %sec

allPt_avgDay = zeros(1440,numel(pt));

for i = 1:numel(pt)
%     i = 12;  %%%%%TEMPORARY
    fs = 400;

    %load in the feature .mat file of interest
    labelLL = 'LL_allCh_2Months_Scaled';
    labelEnergy = 'Energy_allCh_2Months_Scaled';
    labelNLenergy = 'NLenergy_allCh_2Months_Scaled';
    labelArea = 'Area_allCh_2Months';
    labelAlpha = 'alpha_allCh_2Months';
    labelBeta  = 'beta_allCh_2Months';
    labelGamma = 'gamma_allCh_2Months';
    labelTheta = 'theta_allCh_2Months';
    labelHurst = 'hurstExp_allCh_2Months';
    
    allLabels = {labelLL labelEnergy labelNLenergy labelArea labelAlpha ...
        labelBeta labelGamma labelTheta labelHurst};
    
%%%%%%%%%%%%%%%%%%%%%%%% SEIZURE SPIKES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% This section will create a plot of the extracted features with the
% different type of seizure over layed onto
    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate

    
    %load in all features
%     curMat = [pt{i} '_' labelEnergy '.mat'];
%     load(curMat);
    allFeats = cell(length(allLabels),1);
    for fts = 1:length(allLabels)
        curMat = [pt{i} '_' allLabels{fts} '.mat'];
        tmp = cell2mat(struct2cell(load(curMat)));
        %%%
        % Normalize the average feature across channels
        tmp = mean(tmp,2);
        
        avgFeat = mean(tmp);
        stdFeat = std(tmp);
        
        %normalize
        tmp = (tmp - avgFeat)./stdFeat;
        tDays = (1:length(tmp))/60/24;
        maxFeat = mean(tmp) + 8*std(tmp);
        plot(tDays,tmp);
        hold on;
        ylim([0,maxFeat])
        xlabel('Time (days)')
        ylabel('Feature (units vary)')
        
        Annots = [];
    timeUS = [];
    channels = [];
    szTimes = [];
    c = ['k','g','m','y','b'];
    %check to see if current pt. has annotations
    if (size(session.data.annLayer,1) ~= 0)
        for layer = 1:size(session.data.annLayer,2)
           layerName = session.data.annLayer(layer).name;
           [~, timesUS{layer}, ~] = getAllAnnots(session.data,layerName);

           startT = (timesUS{layer}(:,1)*1e-6)/60/60/24;
           startT(startT > 60) = [];

           if startT
                vline(startT,c(layer));
           else
               %don't plot because no startT
           end

            szTimes = [szTimes; startT*24*60];
        end

    else

        %do nothing if there are no annotations
    end

    %overlay the feature values over time centered around the seizure start
    dispWin = 50;
    winSt = szTimes - dispWin;
    winEnd = szTimes + dispWin;
    szT = -dispWin:dispWin;
    figure(20)
    hold on;
    for sz = 1:length(winSt)
        plot(szT,tmp(winSt(sz):winEnd(sz)))
    end
    vline(0)
        
%         allFeats{fts} = cell2mat(struct2cell(load(curMat)));
    end
    
    
   
    

    
    
    
    
    % Now we need to do some kind of automated detection to see if there are
    % similarities between spikes/humps and seizures



%%%%%%%%%%%%%%%%%%%%%%%% SLEEP CYCLE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% This section will create a plot of overlayed features to see if the hump
% trend persists and is a function of sleep patterns

%     curMat = [pt{i} '_' labelNLenergy '.mat'];
%     load(curMat);
% 
%     avgFeat = mean(feat,2);
%     trueMean = mean(avgFeat);
%     avgFeat(avgFeat == 0) = mean(avgFeat);
%     %normalize feature 0 mean unit variance
%     avgFeat = (avgFeat - mean(avgFeat))./std(avgFeat);
%     
%     tDays = (1:length(avgFeat))/60/24;
%     maxFeat = mean(avgFeat) + 2*std(avgFeat);
%     minFeat = mean(avgFeat) - 2*std(avgFeat);
% 
% %     plot(tDays,avgFeat);
% %     hold on;
% %     ylim([0,maxFeat])
% %     xlabel('Time (days)')
% %     ylabel('Feature (units vary)')
% 
%     days = 60;
%     numIdx = length(avgFeat)/days;
% 
%     avgDay = reshape(avgFeat,numIdx,days);
%     avgFeat_Day = mean(avgDay,2);
%     stdFeat_Day = std(avgDay,[],2);
%     timeHours = linspace(0,24,numIdx)';
%     
%     allPt_avgDay(:,i) = avgFeat_Day;
% %     figure(i)
% %     xlim([0,24])
% %     ylim([0,maxFeat])
% %     hold on;
% %     plot(repmat(timeHours,1,60),avgDay,'b','LineWidth',0.5)
% %     plot(timeHours,avgFeat_Day,'k','LineWidth',3);
%     
%     figure(i)
%     errorbar(timeHours,avgFeat_Day,stdFeat_Day);
%     hold on;
%     plot(timeHours,avgFeat_Day,'k','LineWidth',3);
%     ylim([minFeat,maxFeat])
%     xlim([0,24])
%     set(gca,'FontSize',20);
%     set(gca,'LineWidth',2);
%     set(gcf,'Position',get(0,'Screensize')); 
%     xlabel('Time (Hours)')
%     ylabel('Normalized Feature (Non-Linear Energy)')
% %     title([pt{i} ' Normalized One Minute Windowed Non-Linear Energy Averaged Across Days'])
%     
%     
%     labelPlot = [pt{i} 'nlEnergy_dayCycle_noTitle'];
%     print(labelPlot,'-dpng');
%     close;
% 

end

%%
% code is used to average the day feature across all pt. the problem with
% this is that each pts. sleep cylce is different therfore the "dip" occurs
% in a different location each graph would need to be centered in order for
% this graph to be accurate


% avgFeat_allPt = mean(allPt_avgDay,2);
% stdFeat_allPt = std(allPt_avgDay,[],2);
% 
% figure(22)
% errorbar(timeHours,avgFeat_allPt,stdFeat_allPt);
% hold on;
% plot(timeHours,avgFeat_allPt,'k','LineWidth',3);
% ylim([minFeat,maxFeat])
% xlim([0,24])
% set(gca,'FontSize',15);
% set(gca,'LineWidth',2);
% set(gcf,'Position',get(0,'Screensize')); 
% xlabel('Time (Hours)')
% ylabel('Normalized Feature (Non-Linear Energy)')
% title(['All Paitents Normalized Non-Linear Energy Averaged Over Each Day'])
