
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))


session = IEEGSession('NVC1001_23_002','jaredwil','jar_ieeglogin.bin') ;
fs = session.data.sampleRate;               %Find sampling Rate

%Get power and ll in one second windows set in 1 day increments
day = 86400; %sec
hour = 3600; %sec
min = 60; %sec;

%not par samething
tic;
[llN, numNanN] = calcFeature_wil(session.data,1 ,'ll',30,'parTest',[0 10*day], hour,  0);
nPar = toc

%par attempt
tic;
[ll, numNan] = calcFeature_wil(session.data,1 ,'ll',30,'parTest',[0 10*day], hour,  1);
Par = toc


