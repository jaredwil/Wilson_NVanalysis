function [ wk1Feats,wk2Feats ] = diurnalPlot(ptSession)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


    %define all time in seconds
    day = 86400; %sec
    hour = 3600; %sec
    min = 60; %sec

    
    week = day*7;
    weekAust = ceil(week*(400/399.6));

    Norm1  =  [1 + 4*week week + 5*week];
    Norm2  =  [1 +50*week week + 51*week];
    Aust1  =  [1 + 4*weekAust weekAust + 5*weekAust];
    Aust2  =  [1 +50*weekAust weekAust + 51*weekAust];
    
    ch = 1:16;

    [tmp1, numN1] = calcFeature_NV(ptSession.data, ch ,'ll', min, 1 ,labelArea, Norm1 , hour,  0);
    [tmp2, numN2] = calcFeature_NV(ptSession.data, ch ,'ll', min, 1 ,labelArea, Norm2 , hour,  0);
    
    [tmpA1, numA1] = calcFeature_NV(ptSession.data, ch ,'ll', min, 1 ,labelArea, Aust1 , hour,  0);
    [tmpA2, numA1] = calcFeature_NV(ptSession.data, ch ,'ll', min, 1 ,labelArea, Aust2 , hour,  0);


    allTmp = {tmp1 tmp2 tmpA1 tmpA2};
    
    for i = 1:numel(allTmp)
        tmpFeat = allTmp{i};
        
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

