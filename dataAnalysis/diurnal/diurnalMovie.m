%Jared D. Wilson
%7/29/2015
%Create a movie of the diurnal plot for a feature for one channel at 400 Hz
%vs a movie for one channel at 399.6 Hz.

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
dayAust = day*(400/399.6097561);
hour = 3600; %sec
hourAust = 3600*(400/399.6097561);
minute = 60; %sec
minAust = 60*(400/399.6097561);
week = day*7;
weekAust = week*(400/399.6097561);

ch = 1;

numWks = 52;
winLen = minute;
tmpLen = week/minute;
avgFeat_Day = zeros((day/minute),numWks);
stdFeat_Day = zeros((day/minute),numWks);
% feat = zeros((numWks*week*fs/winLen),numel(ch));
%go through an entire year
for wk = 1:numWks
    weekSt  = ((wk-1)*week) + 1;
    weekEnd = (wk)*week;
    label = ['week' num2str(wk)];

    %get the current week
    [tmp, numNan] = calcFeature_NV(session.data, ch ,'ll', minute, 1 ,label, [weekSt weekEnd] , hour,  1); 
%     [~, ~, ~, ~, tmp,numNan] = calcBandPower(session.data, ch ,'gamma',minute,label,[weekSt weekEnd], hour,  1);

%     feat((tmpLen*(wk-1)+1):tmpLen*(wk),:) = tmp;
    
    if(size(tmp,2) > 1)
       tmp = mean(tmp,2); 
    end
    
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

timeHours = linspace(0,24,numIdx)';
%make a movie
for wk = 1:52
    %create graph
    figure(1)
    errorbar(timeHours,avgFeat_Day(:,wk),stdFeat_Day(:,wk));
    plot(timeHours,avgFeat_Day(:,wk))
    xlabel('Time (h)')
    ylim([-2.5 2.5])
    title(['Diurnal Pattern Week: ' num2str(wk)])
    xlim([0 24])
    %wait for a couple seconds
    pause(1)
end

