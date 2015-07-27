%Jared Wilson
%7/27/2015
%
%This script is used with the diurnal plot function to create plots to
%compare the dirunal cylce of features extracted at varrying sampling
%frequencies
%
%
%Start 

clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
%add path to all mat file with feature file in them.
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\allCh_2months_OneminWinFeats'))
set(0, 'DefaulttextInterpreter', 'none') 

pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%define all time in seconds
day = 86400; %sec
hour = 3600; %sec
min = 60; %sec

i = 4;

session = IEEGSession(pt{i},'jaredwil','jar_ieeglogin.bin') ;
fs = session.data.sampleRate;               %Find sampling Rate

