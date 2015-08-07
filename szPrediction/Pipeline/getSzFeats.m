function [extFeats, predIdx] = getSzFeats(ptSession, szStartT, szEndT, szHorizon, winLen, winDisp)
%[trainFeats, predIdx] = getSzFeats(ptSession, szStartT, szEndT, szHorizon, winLen, winDisp)
%   This function takes in the current patient iEEG session and extracts
%   windowed features based on winLen and winDisp entered by the user.
%   Features are only extracted from times that are within the determined
%   szHorizon (entered by user) and the seizure start times (szStartT).
% 
%
% Inputs:
%           ptSession - Entered IEEG session of pt. of interest
%           szStartT  - An array of times (seconds) that note the start of
%           a determined seizure
%           szHorizon - The time before a seizure the is of interest and
%           will be used to predict upcoming sz. (hours)
%           winLen    - Window Length for feature extraction
%           winDisp   - Window Displacement between sucessive windows
%
% Outputs:
%           feats     - calculated feature vector (cell)
%           predIdx   - contains two columns containing the start of the
%           seizure prediction horizon and 30 mins after end of sz. 
%
%%

%%
%anonymous functions
CalcNumWins = @(xLen, fs, winLen, winDisp)floor((xLen-(winLen-winDisp)*fs)/(winDisp*fs));

%%
%Define constant variables
numSz = length(szStartT);  %number of seizures in cur Pt.

numCh = length(ptSession.data.channels);
fs = ptSession.data.sampleRate;


%THIS VALUE CAN CHANGE IS THE NUMBER OF FEATURES EXTRACTED FROM EACH
%CHANNEL
numFeats = size(FeatExt(ptSession.data.getvalues(1:200,1:numCh),fs),2);
%%
%Begin Function
bufferT = 30*60; %buffer time (seconds) after sz considered interictal currently set to 30mins

predIdx = zeros(length(szStartT),2);  %array containing all Idx to be ignored when finding inter-ictal data for training
extFeats = cell(numSz);
%Loop through all Sz.
for sz = 1:numSz
    horizStart = szStartT(sz) - szHorizon*3600; %compute the horizon start time1
    %compute buffered end time of seizure
    if(sz > 1)
        endBufferT = szEndT(sz-1) + bufferT; %add 30 mins in seconds to seizure end time
    else
        endBufferT = 0;
    end
    
    
    %check to see if start of SzHorizon before the end the the sz prior + buffer.
    if(horizStart < endBufferT)
        horizStart = endBufferT;       
    end 
            
    %Get DA DATA! 
    err = 0;
    %try to get data ten times this is to prevent timeouts
    while err < 4
        try    
            dataPreSz = ptSession.data.getvalues((horizStart*fs):(szStartT(sz)*fs),1:numCh);
            break;
        catch
            err = err + 1;    
        end
    end
    
    predIdx(sz,1) = horizStart*fs;    %Keep track of start of sz pred horizon
    predIdx(sz,2) = (szEndT(sz)+bufferT)*fs; %record buffer time as 30 mins after end of sz
    
    dataPreSzNan = dataPreSz; %copy block data with Nans 
    dataPreSz(isnan(dataPreSz)) = 0; %NaN's turned to zero

    numWins = CalcNumWins(size(dataPreSz,1), fs, winLen, winDisp);
    
    feats = zeros(numWins,numFeats);
    for n = 1:numWins
        %off set included
        startWinPt = round(1+(winDisp*(n-1)*fs));
        endWinPt = round(min([winDisp*n*fs,size(dataPreSz,1)]));               
        y = dataPreSz(startWinPt:endWinPt,:);
        
        numNan = sum(isnan(y),1); 

        %only compute the feature vector for windows that have a
        %reasonable amount of good data.
        if(numNan(1) < 400) 
            %compute feature vector for current window
            feats(n,:) = FeatExt(y,fs);
        else
            feats(n,:) = zeros(1,numFeats);
        end

    end

    labelTime = (szStartT(sz)-horizStart-winLen):(-winDisp):0;
    
    %create cell array to return feature vector and lables
    extFeats{sz} = [labelTime' feats];
    disp(['Progress: ' num2str(sz) '/' num2str(numSz)])
end
    %convert to matrix
    extFeats = cell2mat(extFeats);
    %remove all feature vectors that contain all zeros
    extFeats(sum(extFeats(:,2:end),2) == 0, :) = []; 
end