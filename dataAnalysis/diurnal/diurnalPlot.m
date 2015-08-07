function [ wk1Feats,wk2Feats ] = diurnalPlot(session)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    %define all time in seconds
    day = 86400; %sec
    dayAust = day*(399.6097561/400);

    hour = 3600; %sec
    hourAust = 3600*(399.6097561/400);
    minute = 60; %sec
    minAust = 60*(399.6097561/400);
    week = day*7;
    weekAust = week*(399.6097561/400);

    Norm1  =  [(1 + 30*day) (40*day)];
    Norm2  =  [(1 + 300*day) (310*day)];
    Aust1  =  [(1 + 30*dayAust) (41*dayAust)];
    Aust2  =  [(1 + 300*dayAust) (311*dayAust)];
    
    label = 'US';
    label2 = 'Aust';
    ch = 1:16;

%     [tmp1, numN1] = calcFeature_NV(session.data, ch ,'energy', minute, 1 ,label, Norm1 , hour,  1);
%     [tmp2, numN2] = calcFeature_NV(session.data, ch ,'energy', minute, 1 ,label, Norm2 , hour,  1);

%use delta power
    [~, ~, ~, ~, tmp1,numNan1] = calcBandPower(session.data, ch ,'gamma',minute,label,Norm1, hour,  1);
    [~, ~, ~, ~, tmp2,numNan2] = calcBandPower(session.data, ch ,'gamma',minute,label,Norm2, hour,  1);
    
%     [tmpA1, numA1] = calcFeature_NV(session.data, ch ,'energy', minAust, 1 ,label2, Aust1 , hourAust,  1);
    [~,~,~,~,tmpA1, numA1] = calcBandPower(session.data, ch ,'gamma',minAust,label2,Aust1, hourAust,  1);
    tmpA1 = tmpA1(1:length(tmp1),:);
    numA1 = numA1(1:length(tmp1),:);
    
    
%     [tmpA2, numA2] = calcFeature_NV(session.data, ch ,'energy', minAust, 1 ,label2, Aust2 , hourAust,  1);
    [~,~,~,~,tmpA2, numA2] = calcBandPower(session.data, ch ,'gamma',minAust,label2,Aust2, hourAust,  1);
    tmpA2 = tmpA2(1:length(tmp1),:);
    numA2 = numA2(1:length(tmp1),:);
    
    
%     session.data.setResample(399.6)
%     [tmpA1, numA1] = calcFeature_NV(session.data, ch ,'ll', min, 1 ,label2, Norm1 , hour,  1);
%     [tmpA2, numA1] = calcFeature_NV(session.data, ch ,'ll', min, 1 ,label2, Norm2 , hour,  1);
%     session.data.resetResample;
    
%     usTmp   = {tmp1 tmp2};
%     austTmp = {tmpA1 tmpA2};
    allTmp = {tmp1 tmpA1 tmp2 tmpA2};
    subTits = {'Week 3 400Hz'  'Week 3 399.6Hz' 'Week 50 400Hz' 'Week 50 399.6Hz'};
    
    figure(1)
    
    %do computation for sampling rate of 400Hz
    for i = 1:numel(allTmp)
        tmpFeat = allTmp{i};
        
        avgFeat = mean(tmpFeat,2);
        
        %normalize feature 0 mean unit variance
        avgFeat = (avgFeat - mean(avgFeat))./std(avgFeat);
        
        maxFeat = mean(avgFeat) + 2*std(avgFeat);
        minFeat = mean(avgFeat) - 2*std(avgFeat);

        days = ceil(diff(Norm1)/60/60/24);

        numIdx = floor(length(avgFeat)/days);

        
        avgDay = reshape(avgFeat,numIdx,days);
        avgFeat_Day = mean(avgDay,2);
        stdFeat_Day = std(avgDay,[],2);
        timeHours = linspace(0,24,numIdx)';
        
        subplot(2,2,i)
        errorbar(timeHours,avgFeat_Day,stdFeat_Day);
        hold on;
        plot(timeHours,avgFeat_Day)
        xlabel('Time (h)')
        ylim([minFeat maxFeat])
        title(subTits{i})
        xlim([0 24])
    end
   
    
    
end

