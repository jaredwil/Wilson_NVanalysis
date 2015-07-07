%This script is used to find the LL/Energy over the first 2 Months for each
%patient

% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))

pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%define in seconds
day = 86400; %sec
hour = 3600; %sec
min = 60; %sec;

alphaP = cell(length(pt),1);
betaP = cell(length(pt),1);
gammaP = cell(length(pt),1);
thetaP = cell(length(pt),1);

numNan = cell(length(pt),1);
ch = 1:16;

for i = 1:length(pt)
    disp(['Progress: ' num2str(i) '/14'])

    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate

    labelAlpha = 'alphaBP_allCh_2Months';
    labelBeta  = 'betaBP_allCh_2Months';
    labelGamma = 'gammaBP_allCh_2Months';
    labelTheta = 'thetaBP_allCh_2Months';
    
%    [alphaP{i}, numNan{i}] = calcBandPower(session.data, ch ,'alpha',min,labelAlpha,[0 60*day], hour,  1);
%    [betaP{i}, numNan{i}]  = calcBandPower(session.data, ch ,'beta',min,labelBeta,[0 60*day], hour,  1);
     [gammaP{i}, numNan{i}] = calcBandPower(session.data, ch ,'gamma',min,labelGamma,[0 60*day], hour,  1);
    %[thetaP{i}, numNan{i}] = calcBandPower(session.data, ch ,'theta',min,labelTheta,[0 60*day], hour,  1);

end
