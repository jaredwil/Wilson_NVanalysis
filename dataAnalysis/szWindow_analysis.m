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
    
    plotLabels = {'Line Length' 'Energy' 'Non-Linear Energy' 'Area' ...
        'Alpha Power' 'Beta Power' 'Gamma Power' 'Theta Power' ...
        'Hurst Exponent'};
    

%%
% This section will create a plot of the extracted features with the
% different type of seizure over layed onto
    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate

    
    %load in all features
%     curMat = [pt{i} '_' labelEnergy '.mat'];
%     load(curMat);
    allFeats = cell(length(allLabels),1);
%     for fts = 1:length(allLabels)
        fts = 1;
        curMat = [pt{i} '_' allLabels{fts} '.mat'];
        tmp = cell2mat(struct2cell(load(curMat)));
        %%%
        % Normalize the average feature across channels
        tmp = mean(tmp,2);
        avgFeat = sum(tmp)./sum(tmp~=0);
        stdFeat = std(tmp);
        %normalize %only nonzero stuffs
        noZeroIdx = tmp~=0;
        tmp(noZeroIdx) = (tmp(noZeroIdx) - avgFeat)./stdFeat;
        
        if(strfind(allLabels{fts},'hurst')) %don't normalize hurst exp
            %do nothing
        else
        avgFeat = sum(tmp)./sum(tmp~=0);
        stdFeat = std(tmp);
        %normalize %only nonzero stuffs
        noZeroIdx = tmp~=0;
        tmp(noZeroIdx) = (tmp(noZeroIdx) - avgFeat)./stdFeat;
        end
        
        
        tDays = (1:length(tmp))/60/24;
        maxFeat = mean(tmp) + 8*std(tmp);
        figure(1)
        plot(tDays,tmp);
        hold on;
        ylim([0,maxFeat])
        xlabel('Time (days)')
        ylabel('Feature (units vary)')
        
        Annots = [];
        timeUS = [];
        channels = [];
        szTimes = [];
        szendTimes = [];
        c = distinguishable_colors(size(session.data.annLayer,2));
        %check to see if current pt. has annotations
        if (size(session.data.annLayer,1) ~= 0)
            for layer = 1:size(session.data.annLayer,2)
               layerName = session.data.annLayer(layer).name;
               if(isempty(strfind(lower(layerName),'seizure')))
                   continue
               end
               
               [~, timesUS{layer}, ~] = getAllAnnots(session.data,layerName);

               startT = (timesUS{layer}(:,1)*1e-6)/60/60/24;
               endT = (timesUS{layer}(:,2)*1e-6)/60/60/24;
               startT(startT > 60) = [];
               endT(endT > 60) = [];

               if startT
%                   vline(startT,'Color',c(layer,:));
               else
                   %don't plot because no startT
               end

                szTimes = [szTimes; startT*24*60];
                szendTimes = [szendTimes; endT*24*60];
            end

        else

            %do nothing if there are no annotations
        end

%%%%%%%%%%%%%%%%%%%%%%%% SEIZURE SPIKES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% This section gets all of the seziure windows and averages the feature
% across all pt. to view the result. 

        dispWin = 50;
        winSt = szTimes - dispWin;
        winEnd = szTimes + dispWin;
        winTmp = zeros(numel(winSt),dispWin*2+1);
        %overlay the feature values over time centered around the seizure start
        

        szT = -dispWin:dispWin;
        %figure(fts+20)
        %hold on;
        for sz = 1:length(winSt)
%             plot(szT,tmp(winSt(sz):winEnd(sz)))
            winTmp(sz,:) = tmp(winSt(sz):winEnd(sz));
        end
        vline(0)
        
        avgWin(fts,i,:) = sum(winTmp,1)./sum(bsxfun(@plus,bsxfun(@times,winTmp,stdFeat),avgFeat)~=0,1);
        %need to unormalized 
%         if(strfind(allLabels{fts},'hurst')) %don't don't use avg and std to unormalize
%             avgWin(fts,i,:) = sum(winTmp,1)./sum(winTmp~=0,1);
%         else
%             avgWin(fts,i,:) = sum(winTmp,1)./sum(bsxfun(@plus,bsxfun(@times,winTmp,stdFeat),avgFeat)~=0,1);
%         end
%         figure(fts+20)
%         plot(szT,avgWin);
        %allFeats{fts} = cell2mat(struct2cell(load(curMat)));
    end
    
    
   

    %eliminate nans that resulted from dividing by zero
    avgWin(isnan(avgWin)) = 0;
    
    featAvgs_allPt = zeros(length(allLabels),size(avgWin,3));
    featAvgs_allPt = squeeze(sum(avgWin,2)./sum(avgWin~=0,2));
    featStd_allPt = squeeze(std(avgWin,[],2));

    for ft = 1:length(allLabels)
        figure(ft)
        errorbar(szT,featAvgs_allPt(ft,:),featStd_allPt(ft,:))
        hold on;
        plot(szT,featAvgs_allPt(ft,:))
        axis('tight')

        vline(0)
        hline(0,'k')
        set(gca,'FontSize',20);
        set(gca,'LineWidth',2);
        set(gcf,'Position',get(0,'Screensize')); 
        xlabel('Time (Minutes)')
        ylabel(['Normalized ' plotLabels{ft}])
        title(['Normalized ' plotLabels{ft} ' Averaged Across All Pt. for 50 Minutes Surrounding Seizures'],'FontSize',10)

%         labelPlot = [plotLabels{ft} , '_avgAllPt'];
%         print(labelPlot,'-dpng');
%         close;

    end
    
 
   