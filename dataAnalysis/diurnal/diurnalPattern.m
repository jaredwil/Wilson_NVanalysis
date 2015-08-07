%Jared Wilson

%clear the workspace and console
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
%     i = 9;  %%%%%TEMPORARY
    fs = 400;
    
    fsAust = 399.6097561;

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
    
    plotLabels = {'Line Length' 'Energy' 'Non-Linear Energy' 'Area' ...
        'Alpha Power' 'Beta Power' 'Gamma Power' 'Theta Power' ...
        'Hurst Exponent'};
    


%%%%%%%%%%%%%%%%%%%%%%%% SLEEP CYCLE ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
% This section will create a plot of overlayed features to see if the hump
% trend persists and is a function of sleep patterns

    curMat = [pt{i} '_' allLabels{3} '.mat'];
    load(curMat);

    avgFeat = mean(feat,2);
    trueMean = mean(avgFeat);
    avgFeat(avgFeat == 0) = mean(avgFeat);
    %normalize feature 0 mean unit variance
    avgFeat = (avgFeat - mean(avgFeat))./std(avgFeat);
    
    
    
    
%     tDays = (1:length(avgFeat))/60/24;
    
    maxFeat = mean(avgFeat) + 2*std(avgFeat);
    minFeat = mean(avgFeat) - 2*std(avgFeat);

%     plot(tDays,avgFeat);
%     hold on;
%     ylim([0,maxFeat])
%     xlabel('Time (days)')
%     ylabel('Feature (units vary)')

    days = 60;
    numIdx = length(avgFeat)/days;

    avgDay = reshape(avgFeat,numIdx,days);
    avgFeat_Day = mean(avgDay,2);
    stdFeat_Day = std(avgDay,[],2);
    timeHours = linspace(0,24,numIdx)';
    
    allPt_avgDay(:,i) = avgFeat_Day;
%     figure(i)
%     xlim([0,24])
%     ylim([0,maxFeat])
%     hold on;
%     plot(repmat(timeHours,1,60),avgDay,'b','LineWidth',0.5)
%     plot(timeHours,avgFeat_Day,'k','LineWidth',3);
    
    figure(i)
    errorbar(timeHours,avgFeat_Day,stdFeat_Day);
    hold on;
    plot(timeHours,avgFeat_Day,'k','LineWidth',3);
    ylim([minFeat,maxFeat])
    xlim([0,24])
%     set(gca,'FontSize',15);
%     set(gca,'LineWidth',2);
%     set(gcf,'Position',get(0,'Screensize')); 
    xlabel('Time (Hours)')
    ylabel('Normalized Feature (Non-Linear Energy)')
    title([pt{i} ' Normalized One Minute Windowed Non-Linear Energy Averaged Across Days'])
    
    %labelPlot = [pt{i} 'nlEnergy_dayCycle_noTitle'];
%     print(labelPlot,'-dpng');
%     close;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%% AUSTRALIAN STUFFF %%%%%%%%%%%%%%%%%%%%%%%%%%%
    tDaysAust = (1:length(avgFeat))/(1438.6);
    austDayidx = ceil(tDaysAust);
    
    %sort feature matrix into australian days cell array
    austFeats = cell(60,1);
    for d = 1:60
       austFeats{d} = avgFeat(austDayidx == d)';
       if(length(austFeats{d}) > 1438)
          austFeats{d}(end) = []; 
       end
    end
    
    austFeats = cell2mat(austFeats);
    
    avgaustFeat_Day = mean(austFeats,1);
    stdaustFeat_Day = std(austFeats,[],1);
    timeaustHours = linspace(0,24,length(avgaustFeat_Day))';
    
    figure(i+14)
    errorbar(timeaustHours,avgaustFeat_Day,stdaustFeat_Day);
    hold on;
    plot(timeaustHours,avgaustFeat_Day,'k','LineWidth',3);
    ylim([minFeat,maxFeat])
    xlim([0,24])
%     set(gca,'FontSize',15);
%     set(gca,'LineWidth',2);
%     set(gcf,'Position',get(0,'Screensize')); 
    xlabel('Time (Hours)')
    ylabel('Normalized Feature (Non-Linear Energy)')
    title([pt{i} ' Normalized One Minute Windows AUSTRALIA'])
    


end   %this end corrosponds to going throgh all pts.
   
%%
% code is used to average the day feature across all pt. the problem with
% this is that each pts. sleep cylce is different therfore the "dip" occurs
% in a different location each graph would need to be centered in order for
% this graph to be accurate


% avgFeat_allPt = mean(allPt_avgDay,2);
% stdFeat_allPt = std(allPt_avgDay,[],2);
% 
figure(22)
errorbar(timeHours,avgFeat_allPt,stdFeat_allPt);
hold on;
plot(timeHours,avgFeat_allPt,'k','LineWidth',3);
ylim([minFeat,maxFeat])
xlim([0,24])
set(gca,'FontSize',15);
set(gca,'LineWidth',2);
set(gcf,'Position',get(0,'Screensize')); 
xlabel('Time (Hours)')
ylabel('Normalized Feature (Non-Linear Energy)')
title(['All Paitents Normalized Non-Linear Energy Averaged Over Each Day'])