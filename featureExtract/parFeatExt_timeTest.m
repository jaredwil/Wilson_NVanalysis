
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

%number of days to be tested

testN = [1 2 3 4 5 10 15 20 25 30 35 40 45 50 60 70 80 90 100 125 150 200];

timePar    = zeros(length(testN),1);
timeNotPar = zeros(length(testN),1);
check      = zeros(length(testN),1);
parpool(8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1 day %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(testN)
    tic;
    [llN, numNanN] = calcFeature_wil(session.data,1 ,'ll',30,'parTest',[0 testN(i)*day], hour,  0);
    timeNotPar(i) = toc;

    %par attempt
    tic;
    [ll, numNan] = calcFeature_wil(session.data,1 ,'ll',30,'parTest',[0 testN(i)*day], hour,  1);
    timePar(i) = toc;
    
    %ensure that check is 0 meaning there are no erros between the par ll
    %and notpar ll
    check(i) = sum(ll ~= llN);
    
    disp(['Number of Days in Current Test: ' num2str(testN(i))])
end
