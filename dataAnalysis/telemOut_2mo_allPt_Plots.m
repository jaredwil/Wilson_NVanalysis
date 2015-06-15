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

for ptNum = 1:length(pt)
    
    labelLL = [pt{ptNum} '__LL_allCh_2Months.mat'];
    labelEn = [pt{ptNum} '__Energy_allCh_2Months.mat'];
    
    
    ll = load(labelLL);
    ll = ll.feat;
    energy = load(labelEn);
    energy = energy.feat;

    numNanPt = numNan{ptNum};

    telOutage = cell(size(numNanPt,2),1);
    sumOut_sec = cell(size(numNanPt,2),1);
    ll_repOut = ll;
    
    for i = 1:size(numNanPt,2)
        telOutage{i} = find(numNanPt(:,i) > 0);
        sumOut_sec{i} = cumsum(numNanPt(:,i) > 0)*winSize;
        ll_repOut(telOutage{i},i) = 0;
    end

    %plot first channel
    %Find start and end of outages
    out = telOutage{1};
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

    %draw recangles
    figure(1)
    for i = 1:length(outSize)
       rectangle('Position',[startT(i)/5760,-5,outSize(i)/5760,max(ll(:,1))], ...
           'FaceColor',[0.85 0.85 0.85],'EdgeColor','w');
    end
    hold on;
    plot(timeD,ll_repOut(:,1));
    xlabel('Days')
    ylabel('Feature (Line Length)')
    hold on;
    axis([min(timeD) max(timeD) min(ll_repOut(:,1)) 20])

    figure(2)
    plot(timeD, sumOut_sec{1}/86400)
    xlabel('Days')
    ylabel('Total Outage Time (Days)')
end