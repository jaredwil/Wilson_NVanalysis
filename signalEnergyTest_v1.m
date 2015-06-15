%Signal Degradation Optimal Example For AES

% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))


session = IEEGSession('NVC1001_24_001','jaredwil','jar_ieeglogin.bin') ;
fs = session.data.sampleRate;               %Find sampling Rate

%Get power and ll in one second windows set in 1 day increments
day = 86400; %sec
hour = 3600; %sec
min = 60; %sec;

%number of days to be tested

timePar    = zeros(length(testN),1);
timeNotPar = zeros(length(testN),1);
check      = zeros(length(testN),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 day %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [llN, numNanN] = calcFeature_wil(session.data,1:16 ,'ll',30,'parTest',[0 90*day], 4*hour,  0);
    timeNotPar(i) = toc;

