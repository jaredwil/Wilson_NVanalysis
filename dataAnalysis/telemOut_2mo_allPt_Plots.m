%Jared Wilson
%6/15/15

%This script is used to go through all patients and find the time windows
%where dropouts/outages exist. The corrosponding time windows are set to
%zero and the background if set to grey to visuallize the telemetry
%dropouts. Need for this script is the number of NaN (represent dropouts)
%in each time window. 

clear all
close all
clc

addpath(genpath('NVanalysis_data'))
addpath(genpath('Wilson_NVanalysis'))
%this commmand ensure text in tile will not be interpreted as latex
set(0, 'DefaulttextInterpreter', 'none') 
%All patients
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

load('numNan_allCh_allPt_2mo.mat')
timeS = (1:length(ll)).*15;
timeM = timeS./60;
timeH = timeM./60;
timeD = timeH./24;
winSize = 15;
fs = 400;


for ptNum = 1 %:length(pt)
    
    labelLL = [pt{ptNum} '__LL_allCh_2Months.mat'];
    labelEn = [pt{ptNum} '__Energy_allCh_2Months.mat'];
    
    
    ll = load(labelLL);
    ll = ll.feat;
    energy = load(labelEn);
    energy = energy.feat;

    numNanPt = numNan{ptNum};

    telOutage = cell(size(numNanPt,2),1);
    sumOut_sec = cell(size(numNanPt,2),1);
    llOut = ll;
    energyOut = energy;
    maxLL = zeros(1,size(ll,2));
    maxEn = zeros(1,size(ll,2));
    for i = 1:size(numNanPt,2)
        telOutage{i} = find(numNanPt(:,i) > 0);
        sumOut_sec{i} = cumsum(numNanPt(:,i) > 0)*winSize;
        %Set windows with zeros contained in them to zero
        llOut(telOutage{i},i) = 0;
        energyOut(telOutage{i},i) = 0;
        
        avgLL = mean(llOut,1);
        devLL = std(llOut,1);
        avgEn = mean(energyOut,1);
        devEn = std(energyOut,1);
        
        maxLL(i) = avgLL(i) + 5*devLL(i);
        maxEn(i) = avgEn(i) + 5*devEn(i);
        
        %remove outliers
        llOut(llOut(:,i) > maxLL(i),i) = 0;
        energyOut(energyOut(:,i) > maxEn(i),i) = 0;
    end

    %loop through all ch
    for i = 1:size(llOut,2)
        %Find start and end of outages
        out = telOutage{i};
        x = [0; cumsum(diff(out)~=1)];

        numOut = max(x)+1;
        idxStart = ones(numOut,1);
        idxEnd = ones(numOut,1);
        idxEnd(end) = length(x);

        idxEnd(1:end-1) = find(diff(x) > 0);
        idxStart(2:end) = find(diff(x) > 0) + 1;

        startT = out(idxStart);
        endT = out(idxEnd);
        outSize = (endT - startT) + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT LL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1)
        plot(timeD,llOut(:,i));
        hold on;
        for j = 1:length(outSize)
           rectangle('Position',[startT(j)/5760,-5,outSize(j)/5760,max(ll(:,1))], ...
               'FaceColor',[0.8 0.8 0.8],'EdgeColor','w');
        end
        xlabel('Days')
        ylabel('Feature (Line Length)')
        title(['Line Length Over First 60 Days (Channel ' num2str(i) '/Patient ' pt{ptNum} ')'])
        set(gcf,'Color','w');
        axis([min(timeD) max(timeD) min(llOut(:,i)) maxLL(i) + devLL(i)])

        %save as .png and close
        label = ['ll_2mo_ch'  num2str(i) '_' pt{ptNum}];
        print(label,'-dpng');
        close;
        
%%%%%%%%%%%%%%%%%%%%%%%PLOT Energy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(2)
        plot(timeD,energyOut(:,i));
        hold on;
        for j = 1:length(outSize)
           rectangle('Position',[startT(j)/5760,-5,outSize(j)/5760,max(ll(:,1))], ...
               'FaceColor',[0.8 0.8 0.8],'EdgeColor','w');
        end
        xlabel('Days')
        ylabel('Feature (Line Length)')
        title(['Energy Over First 60 Days (Channel ' num2str(i) '/Patient ' pt{ptNum} ')'])
        set(gcf,'Color','w');
        axis([min(timeD) max(timeD) min(energyOut(:,i)) maxEn(i) + devEn(i)])
        %save as .png and close
        label = ['energy_2mo_ch'  num2str(i) '_' pt{ptNum}];
        print(label,'-dpng');
        close;
        
%%%%%%%%%%%%%%%%%%%% Outage Cumulation Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        figure(3)
        plot(timeD, sumOut_sec{i}/86400)
        hold on;
        xlabel('Days')
        ylabel('Total Outage Time (Days)')
        title(['Cumulative Outage Time Over First 60 Days (Patient ' pt{ptNum} ')'])
        set(gcf,'Color','w');
        

    end
     
    %create legend for figure 3 (cumulative plot)
    legend('Channel 1','Channel 2','Channel 3','Channel 4','Channel 5', ...
        'Channel 6','Channel 7','Channel 8','Channel 9','Channel 10', ...
        'Channel 11','Channel 12','Channel 13','Channel 14','Channel 15',...
        'Channel 16')
    %save figure 3 (cumulative plot)
    label = ['CumulativeOutage_' pt{ptNum}];
    print(label,'-dpng');
    close;
end