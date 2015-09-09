%This script is used to find windowed feats over the first 2 Months for
%each patient

% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))

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


for i = 2:length(pt)
% for i = 4
    disp(['Progress: ' num2str(i) '/14'])
    %Start Session
    session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
    fs = session.data.sampleRate;               %Find sampling Rate

    label = 'filt';
    
    labelTest1 = 'spdTest1';
    labelTest2 = 'spdTest2';

     [hwhm, numNan{i}] = calcFeature_NV(session.data, ch ,'hwhm', min, 1,label,[0 100*day], hour,  1);

%    [energy{i}, numNan{i}] = calcFeature_NV(session.data, ch ,'energy', min, labelEnergy,[0 60*day], hour,  1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST TEST TEST TEST TEST TEST %%%%%%%%%%%%%%%%
%      [test, testNan] = calcFeature_NV(session.data, ch ,'hwhm', min, 1,label,[0 1*hour], hour,  0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% OLD STUFF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = profile('info');
% profile off;
% profsave(p,'parProfile_results')
% save('myparProf.mat','p')

%     labelLL = 'LL_allCh_2Months_Scaled';
%     labelEnergy = 'Energy_allCh_2Months_Scaled';
%     labelNLenergy = 'NLenergy_allCh_2Months_Scaled';
%     labelArea = 'Area_allCh_2Months';
