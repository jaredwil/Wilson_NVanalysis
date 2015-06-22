%Jared Wilson
%6/19/2015
% .mat to .csv conversion for NV .mat files

%This script is inteded to be used for .mat files containing time windowed
%features from the NV dataset. 
close all;
clear all
clc;

%MAKE SURE these values are correct
winSize = 15;  %seconds

addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data'))
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

labelLL = [pt{1} '__LL_allCh_2Months.mat'];

llTest = load(labelLL);
llTest = llTest.feat;
sca = size(llTest,1);

llAvgALL = zeros(sca*length(pt),1);
llStdALL = zeros(sca*length(pt),1);

energyAvgALL = zeros(sca*length(pt),1);
energyStdALL = zeros(sca*length(pt),1);

ptID = [];
ptLoc = [];
timeALL = [];
for ptNum = 1:length(pt)
    curPt = pt{ptNum};
    
    labelLL = [pt{ptNum} '__LL_allCh_2Months.mat'];
    labelEn = [pt{ptNum} '__Energy_allCh_2Months.mat'];
    
    ll = load(labelLL);
    ll = ll.feat;
    energy = load(labelEn);
    energy = energy.feat;
    time = linspace(15,size(ll,1)*15,size(ll,1));
    
    
    llAvg     = mean(ll,2);
    llStd     = std(ll,[],2);
    energyAvg = mean(energy,2);
    energyStd = std(energy,[],2);
   
    llAvgALL((ptNum-1)*sca +1:ptNum*sca)     = llAvg;
    llStdALL((ptNum-1)*sca +1:ptNum*sca)     = llStd;
    energyAvgALL((ptNum-1)*sca +1:ptNum*sca) = energyAvg;
    energyStdALL((ptNum-1)*sca +1:ptNum*sca) = energyStd;
    
    ID = repmat(str2num(curPt(12:14)),1,sca);
    Loc = repmat(str2num(curPt(9:10)),1,sca);
    
    ptID = [ptID ID];
    ptLoc = [ptLoc Loc];
    timeALL = [timeALL time];
    
    
    %write the features for each pt. into their own separate csv file
    featPt = [Loc' ID' time' llAvg llStd energyAvg energyStd];
    
    label = ['feats_' pt{ptNum} '.csv']
    csvwrite(label,featPt)
    
end

ptID = ptID';
ptLoc = ptLoc';
timeALL = timeALL';


% featAll = [ptLoc ptID timeALL llAvgALL llStdALL energyAvgALL energyStdALL];
% csvwrite('featsNV.csv',featAll)



