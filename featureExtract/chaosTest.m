%Test the use of chaosExt

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

ch = 1:16;

for i = 4%:numel(pt)
    disp(['Progress: ' num2str(i) '/14'])
    %Start Session
    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate

    label= 'Eigs_allCh_2Months';

    
     [x, numNan{i}] = chaosExt(session.data, ch ,'hfd', min, 1, label,[0 0.5*day], hour,  0);
%    [energy{i}, numNan{i}] = calcFeature_NV(session.data, ch ,'energy', min, labelEnergy,[0 60*day], hour,  1);

end