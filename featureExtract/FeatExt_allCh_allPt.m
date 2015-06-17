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

ll = cell(length(pt),1);
energy = cell(length(pt),1);
numNan = cell(length(pt),1);
ch = 1:16;


for i = 1:length(pt)
    disp(['Progress: ' num2str(i) '/15'])
    %Start Session
    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate

    labelLL = '_LL_allCh_2Months';
    labelEnergy = '_Energy_allCh_2Months';

    [ll{i}, numNan{i}] = calcFeature_wil(session.data, ch ,'ll',15, labelLL,[0 60*day], 2*hour,  1);
    [energy{i}, numNan{i}] = calcFeature_wil(session.data, ch ,'energy',15, labelEnergy,[0 60*day], 2*hour,  1);


end