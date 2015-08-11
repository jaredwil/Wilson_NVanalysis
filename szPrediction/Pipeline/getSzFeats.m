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

%load in filter 
Num = load('NVlowpass.mat');
Num = Num.Num;

%for parprocessing
datasetFN = ptSession.data.snapName;
username = 'jaredwil'; % <<< These values should also be inputs somehow???
pswd = 'jar_ieeglogin.bin';% <<  ^^
%This finds the number of features extracted from the window
numFeats = size(FeatExt(ptSession.data.getvalues(1:200,1:numCh),fs),2);
%%
%Begin Function
bufferT = 30*60; %buffer time (seconds) after sz considered interictal currently set to 30mins

predIdx = cell(length(szStartT),1);  %array containing all Idx to be ignored when finding inter-ictal data for training
extFeats = cell(numSz);

p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)        %if not p will contain all info about current pool;
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
%initialize paralell pool if available
if(poolsize == 0)

    myCluster = parcluster('local');
    numWork = myCluster.NumWorkers;

    if(numWork <= 2)  %processing is done on a laptop so don't do it in parallel
        parpool('local',1)
    elseif(numWork > 8) %limit the number of workers to 6
        parpool('local',8)
    else  %set up a parallel pool with max number of workers available between 3 and 5
        parpool(myCluster)
    end
end

szEndT = [0; szEndT];
% parfor_progress(numSz);

for sz = 1:numSz
%     warning('off');
%     parSession = IEEGSession(datasetFN,username,pswd); %start session on worker
    
    horizStart = szStartT(sz) - szHorizon*3600; %compute the horizon start time1

    %compute buffered end time of prior seizure to ensure they don't
    %overlap
    endBufferT = szEndT(sz) + bufferT; %add 30 mins in seconds to seizure end time
    
    %check to see if start of SzHorizon before the end the the sz prior + buffer.
    if(horizStart < endBufferT)
        horizStart = endBufferT;       
    end 
            
    %Get DA DATA! 
    err = 0;
    %try to get data 4 times this is to prevent timeouts it request spends
    %to much time in que
    while err < 4
        try    
            dataPreSz = ptSession.data.getvalues((horizStart*fs):(szStartT(sz)*fs),1:numCh);
            break;
        catch
            err = err + 1;    
        end
    end
    
%     predIdx{sz,1} = horizStart*fs;    %Keep track of start of sz pred horizon
%     predIdx{sz,2} = (szEndT(sz)+bufferT)*fs; %record buffer time as 30 mins after end of sz
    
    predIdx{sz} = [ horizStart*fs (szEndT(sz)+bufferT)*fs ];
    
    dataPreSzNan = dataPreSz; %copy block data with Nans 
    dataPreSz(isnan(dataPreSz)) = 0; %NaN's turned to zero
    
    %filter the data
    dataPreSz = filtfilt(Num,1,dataPreSz);

    numWins = CalcNumWins(size(dataPreSz,1), fs, winLen, winDisp);
    
%     feats = zeros(numWins,numFeats);
    feats = cell(numWins,1);
    parfor n = 1:numWins
        %off set included
        startWinPt = round(1+(winDisp*(n-1)*fs));
        endWinPt = round(min([winDisp*n*fs,size(dataPreSz,1)]));               
        y = dataPreSz(startWinPt:endWinPt,:);
        
        numNan = sum(isnan(y),1); 

        %only compute the feature vector for windows that have zero
        %Nans
        if(sum(numNan) < 400 && sum(any(y)) ~= 0) 
            %compute feature vector for current window
%             tmpFeats(n,:) = FeatExt(y,fs);
%             N = N + 1; %increase tracked number of valid windows
            feats{n} = FeatExt(y,fs);
        else 
%             tmpFeats(n,:) = zeros(1,numFeats);
            feats{n} = zeros(1,numFeats);

        end

    end
    
    feats = cell2mat(feats);

    labelTime = (szStartT(sz)-horizStart-winLen):(-winDisp):0;
    
    %create cell array to return feature vector and lables
    extFeats{sz} = [labelTime' feats];
    disp(['Progress: ' num2str(sz) '/' num2str(numSz)])
    
end

%convert to matrix
extFeats = cell2mat(extFeats);
predIdx  = cell2mat(predIdx);
%remove all feature vectors that contain all zeros
extFeats(sum(extFeats(:,2:end),2) == 0, :) = []; 

end