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

for i = 4%:numel(pt)
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


    
    %load in all features
%     curMat = [pt{i} '_' labelEnergy '.mat'];
%     load(curMat);
    allFeats = cell(length(allLabels),1);
    for fts = 1:length(allLabels)
%         fts = 1;
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
        
%         if(strfind(allLabels{fts},'hurst')) %don't normalize hurst exp
%             %do nothing
%         else
%         avgFeat = sum(tmp)./sum(tmp~=0);
%         stdFeat = std(tmp);
%         %normalize %only nonzero stuffs
%         noZeroIdx = tmp~=0;
%         tmp(noZeroIdx) = (tmp(noZeroIdx) - avgFeat)./stdFeat;
%         end
        
        
        tDays = (1:length(tmp))/60/24;
        maxFeat = mean(tmp) + 8*std(tmp);
        figure(1)
        plot(tDays,tmp);
        hold on;
%         ylim([0,maxFeat])
        xlabel('Time (days)')
        ylabel('Feature (units vary)')
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% FEATURE SPIKE DETECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%This section should be used to find all the spikes in the feature data
%then figure out how mnay of these spikes corrospond to seizures and
%compute a % of (#spikes that are SZ/#SZ). Also look into maybe clustering
%the spikes and surrounding data points to see if there is a way to detect
%a spike due to sz and spike due to burst or other activity.

        %tmp is normalized to unit std and 0 mean so:
        thresh = 5*std(tmp);
        g = find(tmp > thresh);
        
%         g = find(diff(tmp)> thresh);  %Find where Difference is greater 
%                                           % than the defined threshold value

        g_idx = find(diff(g) ~= 1 & diff(g) > 4);    %Eliminate consecutive indices that 
                              %that corrospond to the same spike
                              
        spikeIdx = [g(g_idx); g(end)];
                              
                              
%         [spikeAmp,spikeIdx] = findpeaks(tmp,'MaxPeakWidth',10,'MinPeakHeight', ...
%             thresh,'SortStr','descend','NPeaks',100,'MinPeakDistance',30);
%         spikeIdx = spikeIdx(diff(spikeIdx) > 5);
%         
                              
% IT IS VERY OBVIOUS THAT SPIKES DO NOT OCCUR B/C of SZ              
%         %find spike indices that occur during a sz it
%         szSpike = cell(numel(szTimes),1);
%         for sz = 1:numel(szTimes)
%             szSpike{sz} = find(spikeIdx > szTimes(sz) & spikeIdx < szendTimes(sz));
%         end
%         szSpike = cell2mat(szSpike);

            dispWin = 2;
            winSt = spikeIdx - dispWin;
            winEnd = spikeIdx + dispWin;
            winEnd(winSt < 1) = [];
            winSt(winSt < 1) = [];
   
            spikeT = -dispWin:dispWin;
            
            spkWin = [];
%             figure(fts+20)
%             hold on;
            for spike = 1:length(winSt)
%                   plot(spikeT,tmp(winSt(spike):winEnd(spike)))
                  spkWin(spike,:) = tmp(winSt(spike):winEnd(spike));
            end
%             vline(0)
            
            
            [coef, score, latent] = pca(spkWin);
            
            pcOne = score(:,1);
            pcTwo = score(:,2);
            pcThree = score(:,3);
            pcFour = score(:,4);
            
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
            pcThree = (pcThree - mean(pcThree))./std(pcThree);
            pcFour = (pcFour - mean(pcFour))./std(pcFour);

            pc = [];
            pc = [pcOne pcTwo pcThree pcFour];

            %Clustering of spikes using k means and top two pc's
            group = kmeans(pc,2,'Distance','sqeuclidean','Replicates',8);

            figure(10)
            gscatter(pcOne,pcTwo,group,'rg')
            legend('Spike Group 1','Spike Group 2','Location','best')
            title('Spike Waveforms Clustered by Top Two Principal Components')
            
            spkOne_avg = mean(spkWin(group == 1,:),1);
            spkTwo_avg = mean(spkWin(group == 2,:),1);


            %Final plot showing clustering results
            figure(11*i)
            plot(spikeT,spkWin(group == 1,:),'r')
            hold on;
            plot(spikeT,spkWin(group == 2,:),'g')
            plot(spikeT,spkOne_avg,'k','LineWidth',4)
            plot(spikeT,spkTwo_avg,'k','LineWidth',4)
            xlabel('Time (mins)')
            ylabel('Feature Normalized')
            title(['Clustered '  plotLabels{fts} ' Spikes in Human NV Data for ' pt{i}])

%                               
%     end
    
    
   
end
