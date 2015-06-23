
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))


session = IEEGSession('NVC1001_24_001','jaredwil','jar_ieeglogin.bin') ;
fs = session.data.sampleRate;               %Find sampling Rate

%Get power and ll in one second windows set in 1 day increments
day = 86400; %sec
hour = 3600; %sec
min = 60; %sec;

%not par samething
% tic;
% [llN, numNanN] = calcFeature_wil(session.data,1 ,'ll',30,'parTest',[0 10*day], hour,  0);
% nPar = toc

%par attempt
tic;
[ll, numNan] = calcFeature_NV(session.data,1 ,'ll',15*min,'parTest',[0 4*day], 2*hour,  0);
Par = toc;


