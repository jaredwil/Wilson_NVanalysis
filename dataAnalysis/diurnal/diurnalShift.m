%Jared D. Wilson
%8/4/2015
%diurnal shift attempt

%Start 
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))

%add path to all mat file with feature file in them.
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\allCh_2months_OneminWinFeats'))
set(0, 'DefaulttextInterpreter', 'none') 

pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

i = 8;

session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin');
fs = session.data.sampleRate;             %Find sampling Rate

%define all time in seconds
day = 86400; %sec
% dayAust = day*(400/399.6097561);
dayAust = day*(399.6097561/400);

hour = 3600; %sec
hourAust = 3600*(399.6097561/400);
minute = 60; %sec
minAust = 60*(399.6097561/400);
week = day*7;
weekAust = week*(399.6097561/400);

ch = 1;
 


numWks = 52;
winLen = minute;
tmpLen = week/minute;
avgFeat_Day = zeros((day/minute),numWks);
stdFeat_Day = zeros((day/minute),numWks);
%go through an entire year
for wk = 1:numWks
    weekSt  = ((wk-1)*week) + 1;
    weekEnd = (wk)*week;
    label = 'week';

    %get the current week
%     [tmp, numNan] = calcFeature_NV(session.data, ch ,'ll', minute, 1 ,label, [weekSt weekEnd] , hour,  1); 
     [~, ~, ~, ~, tmp,numNan] = calcBandPower(session.data, ch ,'gamma',minute,label,[weekSt weekEnd], hour,  1);

%     feat((tmpLen*(wk-1)+1):tmpLen*(wk),:) = tmp;
    
    if(size(tmp,2) > 1)
       tmp = mean(tmp,2); 
    end
    
    feat(((wk-1)*(week/winLen)+1):(wk*(week/winLen)),:) = tmp;
    
    
    tmp = (tmp - mean(tmp))./std(tmp);
    
        
    maxFeat = mean(tmp) + 2*std(tmp);
    minFeat = mean(tmp) - 2*std(tmp);
    
    days = 7;
    numIdx = floor(length(tmp)/days);

    avgDay = reshape(tmp,numIdx,days);
    avgFeat_Day(:,wk) = mean(avgDay,2);
    stdFeat_Day(:,wk) = std(avgDay,[],2);
     
    disp(['Week: ' num2str(wk) ' Done'])

end

%%%%%%%%%%%%%%%%  COMPARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
%Used for scaling
minDif_wk = ((400-399.6097561)*60*60*24*7)/(400*60); %number of samples that should be shifted every week
minDif_day = (400-399.6097561)*60*60*24/(400*60); %number of samples that should be shifted every day

%create a movie showing how to compensate for the difference in sampling
%rate to ensure dirunal pattern does not shift 

%@400Hz
timeHours = linspace(0,24,numIdx)';
figure(1)
for  wk = 1:numWks
    avgFeat_Dayshift(:,wk) = circshift(avgFeat_Day(:,wk),floor(minDif_day*((wk-1)*7)));
    f = plot(timeHours,avgFeat_Dayshift(:,wk));
    grid on;
    xlabel('Time (h)')
    ylim([-2.5 2.5])
    title(['Diurnal Pattern Week: ' num2str(wk)])
    xlim([0 24])
    pause(1);
end




