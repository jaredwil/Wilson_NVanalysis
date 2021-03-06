%Jared D. Wilson
%8/4/2015
%diurnal spectrogram attempt

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


year = [];
feat = zeros((numWks*(week/minute)),numel(ch));
%go through an entire year
for wk = 1:numWks
    weekSt  = ((wk-1)*week) + 1;
    weekEnd = (wk)*week;
    label = ['week'];

    %get the current week
    [tmp, numNan] = calcFeature_NV(session.data, ch ,'nlenergy', minute, 1 ,label, [weekSt weekEnd] , hour,  1); 
%     [~, ~, ~, ~, tmp,numNan] = calcBandPower(session.data, ch ,'gamma',minute,label,[weekSt weekEnd], hour,  1);

%     feat((tmpLen*(wk-1)+1):tmpLen*(wk),:) = tmp;
    
    if(size(tmp,2) > 1)
       tmp = mean(tmp,2); 
    end
    
    feat(((wk-1)*(week/winLen)+1):(wk*(week/winLen)),:) = tmp;
    
    
%     tmp = (tmp - mean(tmp))./std(tmp);
%         
%     maxFeat = mean(tmp) + 2*std(tmp);
%     minFeat = mean(tmp) - 2*std(tmp);
    
    days = 7;
    numIdx = floor(length(tmp)/days);

    avgDay = reshape(tmp,numIdx,days);
    avgFeat_Day(:,wk) = mean(avgDay,2);
    stdFeat_Day(:,wk) = std(avgDay,[],2);
     
    disp(['Week: ' num2str(wk) ' Done'])

end

load('diurnalLP.mat')
load('diurnalLP2.mat')
feat = (feat - mean(feat))./std(feat);
cycFeat = smooth(filter(dirNum2,1,feat),50);
cycFeat = decimate(cycFeat,60);
[s, w, t] = spectrogram(cycFeat,24,[],1028*2*2,24);

y = fft(x);                            % Compute DFT of x
p = unwrap(angle(y));                  % Phase
f = (0:length(y)-1)/length(y)*100;     % Frequency vector
plot(f,p)
xlabel 'Frequency (arb.)'
ylabel 'Phase (rad)'

figure(1)
plot(cycFeat);
ss


