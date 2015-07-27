function [ wk1Feats,wk2Feats ] = diurnalPlot(ptSession)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


    %define all time in seconds
    day = 86400; %sec
    hour = 3600; %sec
    min = 60; %sec
    minAust = 60*(400/399.6);
    
    week = day*7;
    weekAust = ceil(week*(400/399.6));

    Norm1  =  [(1 + 4*week) (week + 4*week)];
    Norm2  =  [(1 + 50*week) (week + 50*week)];
    Aust1  =  [(1 + 4*weekAust) (weekAust + 4*weekAust + min)];
    Aust2  =  [(1 + 50*weekAust) (weekAust + 50*weekAust + min)];
    
    label = 'US';
    label2 = 'Aust';
    ch = 1:16;

    [tmp1, numN1] = calcFeature_NV(ptSession.data, ch ,'ll', min, 1 ,label, Norm1 , hour,  1);
    [tmp2, numN2] = calcFeature_NV(ptSession.data, ch ,'ll', min, 1 ,label, Norm2 , hour,  1);
    
    [tmpA1, numA1] = calcFeature_NV(ptSession.data, ch ,'ll', minAust, 1 ,label2, Aust1 , hour,  0);
    [tmpA2, numA1] = calcFeature_NV(ptSession.data, ch ,'ll', minAust, 1 ,label2, Aust2 , hour,  0);


    usTmp   = {tmp1 tmp2};
    austTmp = {tmpA1 tmpA2};
    
    %do computation for sampling rate of 400Hz
    for i = 1:numel(usTmp)
        tmpFeat = usTmp{i};
        
        avgFeat = mean(tmpFeat,2);
        
        %normalize feature 0 mean unit variance
        avgFeat = (avgFeat - mean(avgFeat))./std(avgFeat);
        
        maxFeat = mean(avgFeat) + 2*std(avgFeat);
        minFeat = mean(avgFeat) - 2*std(avgFeat);

        days = diff(Norm1)/60/24;
        
        
        numIdx = length(avgFeat)/days;

        
        avgDay = reshape(avgFeat,numIdx,days);
        avgFeat_Day = mean(avgDay,2);
        stdFeat_Day = std(avgDay,[],2);
        timeHours = linspace(0,24,numIdx)';

    end
    %do computation for sampling rate of 399.6Hz
    for i = 1:numel(austTmp)
        tmpFeat = austTmp{i};
        
        avgFeat = mean(tmpFeat,2);
        
        %normalize feature 0 mean unit variance
        avgFeat = (avgFeat - mean(avgFeat))./std(avgFeat);
        
        maxFeat = mean(avgFeat) + 2*std(avgFeat);
        minFeat = mean(avgFeat) - 2*std(avgFeat);

        days = 14;
        numIdx = length(avgFeat)/days;

        avgDay = reshape(avgFeat,numIdx,days);
        avgFeat_Day = mean(avgDay,2);
        stdFeat_Day = std(avgDay,[],2);
        timeHours = linspace(0,24,numIdx)';

    end
    
    
end

