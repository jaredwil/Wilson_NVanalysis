function [ interFeats ] = Interict_test( ptSession, szpredIdx, winLen, winDisp, Nlim )
%[ interFeats ] = InterIct_train( ptSession, szpredIdx, szEndT, winLen, winDisp )
%   This function is used to extract interictal data from the NV patient
%   data to be trained on for interictal classification. The times that
%   included sz and szpred Horizon are inputs to be skipped when searching
%   for valid interictal data.
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
%anonymous functions
CalcNumWins = @(xLen, fs, winLen, winDisp)floor((xLen-(winLen-winDisp)*fs)/(winDisp*fs));

%%
%Define constant variables 
datasetFN = ptSession.data.snapName;

numCh = length(ptSession.data.channels);
fs = ptSession.data.sampleRate;
blockLenSecs = 60*60; %block length to be extracted at a time (in seconds)
numWins = CalcNumWins(blockLenSecs*fs,fs,winLen,winDisp); %number of windows per block

N = 0; %number of valid interictal feature vectors extracted

%THIS VALUE CAN CHANGE IS THE NUMBER OF FEATURES EXTRACTED FROM EACH
%CHANNEL
numFeats = size(FeatExt(ptSession.data.getvalues(1:200,1:numCh),fs),2);
%%
%Begin Function
endSearch = ptSession.data.channels(1).get_tsdetails.getDuration/1e6*fs;  %if you reach this index stop search end of training data;
startSearch = szpredIdx(1,2); %start search at the end of the first testing seizure

numBlocks = ceil( ((endSearch/fs)-(startSearch/fs))/blockLenSecs);  %number of blocks to go through
interFeats = cell(numBlocks,1);

randBlocks = randperm(numBlocks);

%initialize startIdx and endIx which tell what should be extracted
% for j = 1:numBlocks
j = 1;
while(j < numBlocks)
    randj = randBlocks(j); %get a random block 
    
    %get data for current block                (include test set offset)
    startBlockPt = 1+(blockLenSecs*(randj-1)*fs)      + startSearch;
    endBlockPt = min(blockLenSecs*randj*fs,endSearch) + startSearch;
    
    %just skip block if either the start or end block pt is inside horizon
    if(sum(startBlockPt > szpredIdx(:,1)) > 0 && sum(startBlockPt < szpredIdx(:,2)) > 0 || ...
            sum(endBlockPt > szpredIdx(:,1)) > 0 && sum(endBlockPt < szpredIdx(:,1)) > 0)
        %skp this block and make tmpFeats all zeros
        numWins = CalcNumWins(blockLenSecs*fs,fs,winLen,winDisp);
        tmpFeats = zeros(numWins,numFeats); %initialize 
    else %do feature extraction
        blockData = ptSession.data.getvalues(startBlockPt:endBlockPt,1:numCh);
        
        if(sum(isnan(blockData),1) == blockLenSecs*fs)
%             j = j + 24; %skip ahead 24 hours
            tmpFeats = zeros(numWins,numFeats);

        else
            blockNan = blockData; %copy block data with Nans 
            blockData(isnan(blockData)) = 0; %NaN's turned to zero

            tmpFeats = zeros(numWins,numFeats);
            for n = 1:numWins
                %off set included
                startWinPt = round(1+(winDisp*(n-1)*fs));
                endWinPt = round(min([winDisp*n*fs,size(blockData,1)]));               
                y = blockData(startWinPt:endWinPt,:);
                %find the number of Nans in the current window
                numNan = sum(isnan(blockNan(startWinPt:endWinPt,:)),1); 

                %only compute the feature vector for windows that have zero
                %Nans
                if(numNan < 400)
                    %compute feature vector for current window
                    tmpFeats(n,:) = FeatExt(y,fs);
                    N = N + 1; %increase tracked number of valid windows
                else 
                    tmpFeats(n,:) = zeros(1,numFeats);
                end
            end
        end
         
    end
    
    interFeats{j} = tmpFeats;
    
    %STOP if one of two conditions are met:
    %       1.) 5000 total interictal feature vectors have been extracted
    %       2.) end of search index has been reached
    if(N > Nlim)
        break;
    end
    
    disp(['Progress: N = ' num2str(N) ' j = ' num2str(j)])
    j = j + 1;
end
    %convert to matrix
    interFeats = cell2mat(interFeats);
    %remove all feature vectors that contain all zeros
    interFeats(sum(interFeats,2) == 0, :) = [];
end

