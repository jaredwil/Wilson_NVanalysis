%This script is used to find the LL/Energy over the first 2 Months for each
%patient

% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))

profile on

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

%CHANGED TO FIND PROFILE
%profile on;
%for i = 1  %:length(pt)
for i = 4
    disp(['Progress: ' num2str(i) '/14'])
    %Start Session
    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate

    labelLL = 'LL_allCh_2Months_Scaled';
    labelEnergy = 'Energy_allCh_2Months_Scaled';

     [ll{i}, numNan{i}] = calcFeature_NV(session.data, ch ,'ll', min, labelLL,[0 3*day], hour,  1);
%    [energy{i}, numNan{i}] = calcFeature_NV(session.data, ch ,'energy', min, labelEnergy,[0 60*day], hour,  1);

end

%p = profile('info');
%profile off;
%profsave(p,'parProfile_results')
%save('myparProf.mat','p')