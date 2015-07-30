%Jared Wilson
%7/29/2015
%Extract EEG from time window features that exceed threshold

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

for i = 8%:numel(pt)
    
    
    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate

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
%start analysis 

%     for fts = 1:length(allLabels)  %get all feats
        fts = 1;  %get LL
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
        
        %find all spea
        %tmp is normalized to unit std and 0 mean so:
        thresh = 5*std(tmp);
        tmpIdx = find(tmp > thresh);  %get all windows that exceed threshold

        tDays = (1:length(tmp))/60/24;
        maxFeat = mean(tmp) + 8*std(tmp);
        figure(1)
        plot(tDays,tmp);
        hold on;
%         ylim([0,maxFeat])
        hline(thresh)
        xlabel('Time (days)')
        ylabel('Feature (units vary)')
        
        %pick out 20 random spikes that don't come from the same 20 min
        %window
        g_idx = find(diff(tmpIdx) ~= 1 & diff(tmpIdx) > 20);    %Eliminate consecutive indices that 
        spikeIdx = tmpIdx(g_idx);
        sampSpike = datasample(spikeIdx,20);
        
        %get 20 windows that are below the threshold;
        belowIdx = 1:length(tmp);
        belowIdx(spikeIdx) = [];
        sampNotSpike = datasample(belowIdx, 20);
        
        %got to portal and get the corrosponding spike window
        winLen = 60; %seconds
        spikeEEG = zeros(winLen*fs,length(sampSpike));
        nospikeEEG = zeros(winLen*fs,length(sampNotSpike));
        timeMin = linspace(1,60,winLen*fs);
%         figure(2)
%         subplot(2,1,1)
        figure(3)
        title('EEG from 20 60 Second Time Windows with LL above 5 std from the Mean')
        hold on;
        %get the spike windows
        for j = 1:length(sampSpike)
            startT = (winLen*(sampSpike(j) - 1)*fs + 1);  
            endT   = winLen*sampSpike(j)*fs;
            eegTmp = session.data.getvalues(startT:endT,1:16);
            
%             eegTmp = session.data.getvalues(startT:endT,1);

            spikeEEG(:,j) = mean(eegTmp,2);
            plot(timeMin,spikeEEG(:,j) + j*20)
            disp(['Done: ' num2str(j/length(sampSpike))])
        end
        figure(4)
%         subplot(2,1,2)
        title('EEG from 20 60 Second Time Windows with LL below 5 std from the Mean')
        hold on;
        %get the non spike windows
        for j = 1:length(sampNotSpike)
            startT = (winLen*(sampNotSpike(j) - 1)*fs + 1);  
            endT   = winLen*sampNotSpike(j)*fs;
            eegTmp = session.data.getvalues(startT:endT,1:16);
            
%             eegTmp = session.data.getvalues(startT:endT,1);

            nospikeEEG(:,j) = mean(eegTmp,2);
            plot(timeMin,nospikeEEG(:,j) + j*2);
            disp(['Done: ' num2str(j/length(sampNotSpike))])
        end
        
%      end    
    
   
end
