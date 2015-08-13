function [ interFeats ] = train_interPar( ptSession, szpredIdx, winLen, winDisp , Nlim)
%[ interFeats ] = train_interPar( ptSession, szpredIdx, szEndT, winLen, winDisp )
% Parallelization attempt for getting interictal data
%   This function is used to extract interictal data from the NV patient
%   data to be trained on for interictal classification. The times that
%   included sz and szpred Horizon are inputs to be skipped when searching
%   for valid interictal data.
%
%       This function utilizes permutations to select a random subset of
%       interictal data to train on vs the first few blocks avaialable.
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
%anonymous functions
CalcNumWins = @(xLen, fs, winLen, winDisp)floor((xLen-(winLen-winDisp)*fs)/(winDisp*fs));

%%
%Define constant variables 
datasetFN = ptSession.data.snapName;
username = 'jaredwil'; % <<< These values should also be inputs somehow???
pswd = 'jar_ieeglogin.bin';% <<  maybe read from .txt file that would be cool

numCh = length(ptSession.data.channels);
fs = ptSession.data.sampleRate;
blockLenSecs = 60*60; %block length to be extracted at a time (in seconds)
numWins = CalcNumWins(blockLenSecs*fs,fs,winLen,winDisp); %number of windows per block

N = 0; %number of valid interictal feature vectors extracted


%THIS VALUE CAN CHANGE IS THE NUMBER OF FEATURES EXTRACTED FROM EACH
%CHANNEL
tstClip = ptSession.data.getvalues(1:200,1:numCh);
tstClip(isnan(tstClip)) = 0; %NaN's turned to zero
stClip = 1;
while(sum(sum(isnan(tstClip),1)) > 400 && sum(any(tstClip)) == 0) 
    tstClip = ptSession.data.getvalues((stClip*200+1):(stClip+1)*200,1:numCh);
    stClip = stClip + 1;
end

numFeats = size(FeatExt(tstClip,fs),2);
%Begin Function
endSearch = szpredIdx(2,end);  %if you reach this index stop search end of training data;

numBlocks = ceil(endSearch/fs/blockLenSecs);  %nubmer of blocks to go through
interFeats = cell(numBlocks,1);

randBlocks = randperm(numBlocks); %create a random order of the available

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Establish parpool
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

    if(numWork < 2)  %processing is done on a laptop so don't do it in parallel
        parpool('local',1)
    elseif(numWork > 8) %limit the number of workers to 6
        parpool('local',8)
    else  %set up a parallel pool with max number of workers available between 2 and 6
        parpool(myCluster)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for j = 1:numBlocks
j = 0;
while(j < numBlocks)
    totN = cell(poolsize,1);
    parfor p = 1:poolsize   %use all available workers in the pool    
        %establish par session on worker
        Npar = 0;
        parSession = IEEGSession(datasetFN,username,pswd); %start session on worker
        
        randj = randBlocks(j+p); %get a random block 

        %get start and end time for current block
        startBlockPt = 1+(blockLenSecs*(randj-1)*fs);
        endBlockPt = min(blockLenSecs*randj*fs,endSearch);

        %just skip block if either the start or end block pt is inside horizon
        if(sum(startBlockPt > szpredIdx(:,1)) > 0 && sum(startBlockPt < szpredIdx(:,2)) > 0 || ...
                sum(endBlockPt > szpredIdx(:,1)) > 0 && sum(endBlockPt < szpredIdx(:,1)) > 0)

            %skp this block and fill tmpFeats with all zeros (removed later)
            tmpFeats = zeros(numWins,numFeats); 

        else %do feature extraction
            %Get DA DATA! 
            err = 0;
            %try to get data 4 times this is to prevent timeouts it request spends
            %to much time in que
            while err < 4
                try    
                    blockData = parSession.data.getvalues(startBlockPt:endBlockPt,1:numCh);
                    break;
                catch
                    err = err + 1;    
                end
            end

            %if the whole block is empty then skip immediatly
            if(sum(isnan(blockData),1) == blockLenSecs*fs) 
                tmpFeats = zeros(numWins,numFeats);
            else
                blockNan = blockData; %copy block data with Nans 
                blockData(isnan(blockData)) = 0; %NaN's turned to zero

                tmpFeats = zeros(numWins,numFeats);
    %             tmpFeats = cell(numWins,1);
                for n = 1:numWins
                    %off set included
                    startWinPt = round(1+(winDisp*(n-1)*fs));
                    endWinPt = round(min([winDisp*n*fs,size(blockData,1)]));               
                    y = blockData(startWinPt:endWinPt,:);
                    %find the number of Nans in the current window
                    numNan = sum(isnan(blockNan(startWinPt:endWinPt,:)),1); 

                    %only compute the feature vector for windows that have zero
                    %Nans
                    if(sum(numNan) < 400 && sum(any(y)) ~= 0) 
                        %compute feature vector for current window
                        tmpFeats(n,:) = FeatExt(y,fs);
                        Npar = Npar + 1; %increase tracked number of valid windows in paralell loop
                    else 
                        tmpFeats(n,:) = zeros(1,numFeats);

                    end
                end
            end

        end

        interFeats{j+p} = tmpFeats;
        totN{p} = Npar;
    
    end
    %calc total N so far
    N = N + sum(cell2mat(totN));
    %show progress
    j = j + poolsize;
    disp(['Progress: N = ' num2str(N) ' j = ' num2str(j)])
    %STOP if one of two conditions are met:
    %       1.) 5000 total interictal feature vectors have been extracted
    %       2.) end of search index has been reached
    if(N > Nlim)
        break;
    end

end
    disp('Interictal Data Complete!!!')
    %convert to matrix
    interFeats = cell2mat(interFeats);
    %remove all feature vectors that contain all zeros
    interFeats(sum(interFeats,2) == 0, :) = [];
end

