%Jared Wilson
%6/19/2015
%This script is for analysis of NV dataset consisting of 14 pt. This script
%will look at the first 60 days of each pt and find the avg and std of
%features across all 16 channels to observe notisable trends. 

%Correlation between channels could also prove to be interesting. 


%Really only the fourth pt. 'NVC1001_25_005' provided a good plot for the
%energy..... other ones are missing data in the first few days which makes
%it difficult to observe a trend. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;

set(0, 'DefaulttextInterpreter', 'none') 
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data'))

pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

load('numNan_allCh_allPt_2mo.mat')
timeS = (1:length(numNan{1})).*15;
timeM = timeS./60;
timeH = timeM./60;
timeD = timeH./24;
winSize = 15;  %window size was 15 seconds
fs = 400;
numCh = 16;

plotNum = 60; %Number of points for avg plot (e.g. 60 - avg across day, 
              % 120 avg. across half day, ...etc)

ptRange = [4 5 6 9 10 11];

for ptNum = ptRange %1:length(pt)
    
    labelLL = [pt{ptNum} '__LL_allCh_2Months.mat'];
    labelEn = [pt{ptNum} '__Energy_allCh_2Months.mat'];
    
    ll = load(labelLL);
    ll = ll.feat;
    energy = load(labelEn);
    energy = energy.feat;
    numNanPt = numNan{ptNum};  %number Nan for current pt.

    %init. variables used to store averaged features
    telOutage = cell(size(numNanPt,2),1);
    llOut = ll;
    dayLL = zeros(size(ll,2),plotNum);
    dayEnergy = zeros(size(ll,2),plotNum);
    energyOut = energy;
    for i = 1:size(numNanPt,2)
        telOutage{i} = find(numNanPt(:,i) > 0);  %find all windows with Nans
        %Set windows with zeros contained in them to zero ch. by ch.
        tmpLL = ll(:,i);
        tmpLL(telOutage{i}) = 0;
        tmpEnergy = energy(:,i); 
        tmpEnergy(telOutage{i}) = 0;
        
        %normalize        
        tmpLLNa = tmpLL;
        tmpLLNa(tmpLL == 0) = [];
        normAv = mean(tmpLLNa);
        normStd = std(tmpLLNa);

        tmpLL = (tmpLL - normAv)./normStd;
        
        tmpEnergyNa = tmpEnergy;
        tmpEnergyNa(tmpEnergy == 0) = [];
        normAv = mean(tmpEnergyNa);
        normStd = std(tmpEnergyNa);

        tmpEnergy = (tmpEnergy - normAv)./normStd;


        %reshape each channel so that each row represents 1 day
        tmpLL = reshape(tmpLL, length(tmpLL)/plotNum, plotNum);
        tmpEnergy = reshape(tmpEnergy, length(tmpEnergy)/plotNum, plotNum);
        
        for d = 1:size(tmpLL,2)
            tmpDayLL = tmpLL(:,d);
            tmpDayLL = tmpDayLL(tmpDayLL~=0);
            tmpDayEnergy = tmpEnergy(:,d);
            tmpDayEnergy = tmpDayEnergy(tmpDayEnergy~=0);            
            if tmpDayLL
                dayLL(i,d) = mean(tmpDayLL);
                dayEnergy(i,d) = mean(tmpDayEnergy);
            end
        end
    end
    
    %compute average across channels for each day
    ccDayLL     = mean(dayLL,1);
    ccDayEnergy = mean(dayEnergy,1);
    
    devDayLL     = std(dayLL);
    devDayEnergy = std(dayEnergy);
    
    figure(1*ptNum)
    errorbar(ccDayLL, devDayLL)
    set(gca,'FontSize',15);
    set(gca,'LineWidth',2);
    set(gcf,'Position',get(0,'Screensize')); 
    xlabel('Day')
    ylabel('Line Length')
    title(['Average Time Windowed Line Length Across Channels for ' pt{ptNum}])
    xlim([0 60])

    figure(2*ptNum)
    errorbar(ccDayEnergy, devDayEnergy)
    set(gca,'FontSize',15);
    set(gca,'LineWidth',2);
    set(gcf,'Position',get(0,'Screensize')); 
    xlabel('Day')
    ylabel('Channel Normalized Energy')
    title(['Average Time Windowed Energy Across Normalized Channel Features for ' pt{ptNum}])
    xlim([0 60])

end



        
     
    