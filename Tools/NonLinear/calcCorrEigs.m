function [eigVals, numNan] = calcCorrEigs(datasets,channels,winLen,outLabel,runIndex,blockLenSecs,parFlag,varargin)
%[feat, numNan] = calcFeature_wil(datasets,channels,feature,winLen,outLabel,runIndex,blockLenSecs,parFlag,varargin)
%       This function is used to analyze how the eigenvalues of the channel
%       correlation matrix changes over long data sets; able to select
%       parflag to run processes in paralell.
%
% Author: Jared Wilson 
% 4/8/15 v1 working
%
% NOTE: This function is a derivative of calcFeature function developed by
% Hoameng Ung
%Usage: calcFeature(datasets,channels,feature,winLen,outLabel,filtFlag,varargin)
%This function will divide IEEGDataset channels into blocks of
%and within these blocks further divide into winLen. Correlation EigenVals
%will be calculated for each winLen and saved in a .mat matrix.

% Parallel Processing Usage: if parFlag set multiple sessions are
% established to pull data blocks simultaneosuly
%
%

%Input:
%   'dataset'       :   IEEG dataset 
%   'channels'      :   vector of channels
%   'winLen'        :   winLen = vector of windowlengths (s) to calculate
%                       features over
%   'outlabel'      :   Suffix to save features to
%   (datasetname_outlabel.mat)
%   'runIndex'      :   Set the length of analysis [beginSec endSec] in
%   seconds or 'all' will run on entire dataset
%   'blockLenSec'   :   Length of blocks to be pulled from portal at a time
%   can be set to improve efficieny. (seconds)
%   'parFlag'       :   [0/1] 1: parallel processing (parfor) will be used
%   to increase speed

%Output:
%   'eigVals'       :   Calculated Correlation Eigen values for each window
%   'numNan'        :   The number of NaN that is contained in each window
%   in each channel same size as feat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              TO DO LIST 
% - add numPar input to designat number of processor to use
% - try to condense inputs to function (there are too many)
% - do error checking 
% - handel inputs to accout for optional inputs being omited
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%   ERROR CHECKS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Error checks before function

%make sure runIndex given is valid
    %value must be size 2 (start and end) and end must be greater than
    %start OR value can be 'all'

%make sure blockLenSecs is valid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deal with optional inputs 
% ???????????????????????? learn this 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Begin Function

%% Anonymous functions

%% Initialization
%blockLenSecs = 3600; %get data in blocks (use as defult)

%% Parallel Processing
if parFlag    %if flag is set do processing in parallel
    
    for i = 1:numel(datasets)

        %if all selected analyze entire dataset
        if (strcmpi('all',runIndex) == 1) 
            datasetFN = datasets(i).snapName;
            duration = datasets(i).channels(1).get_tsdetails.getDuration;
            fs = datasets(i).channels(1).sampleRate;
            numPoints = duration/1e6*fs;
            %calculate number of blocks
            %numBlocks = CalcNumWins(numPoints,fs,blockLenSecs,blockLenSecs);
            numBlocks = ceil(numPoints/fs/blockLenSecs);
            startSec = 0;
        %Else only run analysis on selected data
        else
            datasetFN = datasets(i).snapName;
            startSec = runIndex(1);
            endSec = runIndex(2);
            duration = endSec - startSec; %endSec minus beginSec
            fs = datasets(i).channels(1).sampleRate;
            numPoints = duration*fs;
            %calculate number of blocks
            %numBlocks = CalcNumWins(numPoints,fs,blockLenSecs,blockLenSecs);
            numBlocks = ceil(numPoints/fs/blockLenSecs);    
        end
        
        %% Feature extraction loop
        %paralell variables
        %Number of workers is determined by the pool that is set. This
        %determines the number of tasks that will be allocated to each
        %worker
        numPar = 6;     
        username = 'jaredwil'; % <<< These values should also be inputs somehow???
        pswd = 'jar_ieeglogin.bin';% <<  ^^

        %determine number of blocks per parallel session
        parBlocks = floor(numBlocks/numPar);
        
        %off set index for each processor
        parOffset = parBlocks*blockLenSecs*fs;
        parFeat = cell(numPar,1);
        parNan = cell(numPar,1);
        
        %%%REC CHANGE
        featTot = cell(parBlocks*numPar,1);
        numNanTot = cell(parBlocks*numPar,1);
        
        parfor p = 1:numPar
                %warning('off')
                %Initilize variable to be used specifically in each proc.
                parfeat = cell(parBlocks,1);
                parnumNan = cell(parBlocks,1);
                
                session = IEEGSession(datasetFN,username,pswd);
                
                %create data variable
                parData = session.data;
        
            for j = 1:parBlocks
                %Get data
                startBlockPt = 1+(blockLenSecs*(j-1)*fs)+ startSec*fs;
                startBlockPt = startBlockPt + parOffset*(p-1);
                endBlockPt = min(blockLenSecs*j*fs,numPoints)+ startSec*fs;
                endBlockPt = endBlockPt + parOffset*(p-1);
                
                err = 0;
                %try to get data ten times this is to prevent timeouts
                while err < 10
                    try    
                        blockData = parData.getvalues(startBlockPt:endBlockPt,channels);
                        break;
                    catch
                        err = err + 1    
                    end
                end

                %%%%%DO a check to see if the entire block is Nan if this is
                %%%%%true then skip this block. Ideally returning an indication
                %%%%%that this row was skipped (return 0's) to speed up
                %%%%%computation. 
                if(sum(isnan(blockData),1) == blockLenSecs*fs) % all 0's
                    nChan = numel(channels);
                    %calculate feature every winLen secs
                    numWins = ceil(size(blockData,1)/(winLen*fs));
                    tmpEig = zeros(numWins,nChan);
                    tmpNan  = ones(numWins,nChan)*(winLen*fs);
                else       %do normal feature extraction
                    blockNan = blockData; %copy block data with Nans 
                    blockData(isnan(blockData)) = 0; %NaN's turned to zero
                    nChan = numel(channels);
                    %calculate feature every winLen secs
                    numWins = ceil(size(blockData,1)/(winLen*fs));
                    tmpEig = zeros(numWins,nChan);
                    tmpNan  = zeros(numWins,nChan);
                    for n = 1:numWins
                        %off set included
                        startWinPt = round(1+(winLen*(n-1)*fs));
                        endWinPt = round(min(winLen*n*fs,size(blockData,1)));
                        %Find number of Nan
                        tmpNan(n,:) = sum(isnan(blockNan(startWinPt:endWinPt,:)),1);

                        idxN = find(isnan(blockNan(startWinPt:endWinPt,1)) == 0);  %find all indices without Nan
                        x = [0; cumsum(diff(idxN)~=1)];                 %find consecutive non Nan values

                        
                        %check to see if there are any outages if there are
                        %not then proceed normally else do computation to
                        %componsate for Nans
                        if max(x) == 0
                            y = blockData(startWinPt:endWinPt,1:nChan);
                            
                            %normalize windows zero mean/unit std dev
                            
                            R = corrcoef(y);
                            tmpEig(n,:) = eig(R)';

       
                            
                        else
                            %number of consecutive outs
                       	    numOut = max(x)+1;                                  
                            idxStart = ones(numOut,1);    %start ind                
                            idxEnd = ones(numOut,1);      %end ind
                            idxEnd(end) = length(x);      %last end point is end

                            %find all start and end index of outage
                            idxEnd(1:end-1) = find(diff(x) > 0);
                            idxStart(2:end) = find(diff(x) > 0) + 1;

                            %Start/end/and length of outage
                            startT = idxN(idxStart);   %values are not really time
                            endT = idxN(idxEnd);
                            outLen = (endT - startT);
                        
                            k = 1;
                            winEig = [];
                            for goodWin = 1:numOut
                                if outLen(goodWin) > 5                
                                
                                tmp = blockData(startT(goodWin):endT(goodWin),1:nChan);
                                R = corrcoef(tmp);
                                winEig(k,:) = eig(R)';
                                k = k + 1;
                                end
                            end
                            tmpEig(n,:) = mean(winEig,1);
                            
                        end
                    end

                end

                parfeat{j} = tmpEig;
                parnumNan{j} = tmpNan;
                %This doesn't really work for par processes
%                 percentDone = 100 * j / numBlocks;
%                 msg = sprintf('Percent done: %3.1f',percentDone); %Don't forget this semicolon
%                 fprintf([reverseStr, msg])
%                 reverseStr = repmat(sprintf('\b'), 1, length(msg));
            end
            
        parFeat{p} = parfeat; %Store the result from isolated process in 
        parNan{p} = parnumNan; %its appropriate location in cell
        end
    %put all cells back together in matrix
      eigVals = [];
      numNan = [];
%     feat = cell(numBlocks,1);
%     numNan = cell(numBlocks,1);
    for tm = 1:numPar
%        %%%%THIS IS BAD CODE!!!! CHANGE THIS!
%        tmpF = parFeat{tm};
%        tmpN = parNan{tm};
%        feat = [feat; tmpF];
%        numNan = [numNan; tmpN];
       featTot(parBlocks*(tm-1)+1 : parBlocks*tm) = parFeat{tm};
       numNanTot(parBlocks*(tm-1)+1 : parBlocks*tm) = parNan{tm}; 
       
    end
    
    fprintf('\n');
    eigVals = cell2mat(featTot);
    numNan = cell2mat(numNanTot);
    save([datasetFN '_' outLabel '.mat'],'eigVals','-v7.3');
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Normal Processing
else          %do normal processing if flag is not set
    for i = 1:numel(datasets)
        %if all selected analyze entire dataset
        if (strcmpi('all',runIndex) == 1) 
            datasetFN = datasets(i).snapName;
            duration = datasets(i).channels(1).get_tsdetails.getDuration;
            fs = datasets(i).channels(1).sampleRate;
            numPoints = duration/1e6*fs;
            %calculate number of blocks
            %numBlocks = CalcNumWins(numPoints,fs,blockLenSecs,blockLenSecs);
            numBlocks = ceil(numPoints/fs/blockLenSecs);
            startSec = 0;
        %Else only run analysis on selected data
        else
            datasetFN = datasets(i).snapName;
            startSec = runIndex(1);
            endSec = runIndex(2);
            duration = endSec - startSec; %endSec minus beginSec
            fs = datasets(i).channels(1).sampleRate;
            numPoints = duration*fs;
            %calculate number of blocks
            %numBlocks = CalcNumWins(numPoints,fs,blockLenSecs,blockLenSecs);
            numBlocks = ceil(numPoints/fs/blockLenSecs);    
        end
        
        %% Feature extraction loop
        eigVals = cell(numBlocks,1);
        numNan = cell(numBlocks,1);
        reverseStr = '';
        for j = 1:numBlocks
            %Get data
            startBlockPt = 1+(blockLenSecs*(j-1)*fs)+ startSec*fs;
            endBlockPt = min(blockLenSecs*j*fs,numPoints)+ startSec*fs;
            blockData = datasets(i).getvalues(startBlockPt:endBlockPt,channels);
            

            
            %%%%%DO a check to see if the entire block is Nan if this is
            %%%%%true then skip this block. Ideally returning an indication
            %%%%%that this row was skipped (return 0's) to speed up
            %%%%%computation. 
            if(sum(isnan(blockData),1) == blockLenSecs*fs) % all 0's
                nChan = numel(channels);
                %calculate feature every winLen secs
                numWins = ceil(size(blockData,1)/(winLen*fs));
                tmpEig = zeros(numWins,nChan);
                tmpNan  = ones(numWins,nChan)*(winLen*fs);
            else       %do normal feature extraction
                blockNan = blockData; %copy block data with Nans 
                blockData(isnan(blockData)) = 0; %NaN's turned to zero
                nChan = numel(channels);
                %calculate feature every winLen secs
                numWins = ceil(size(blockData,1)/(winLen*fs));
                tmpEig = zeros(numWins,nChan);
                tmpNan  = zeros(numWins,nChan);
                for n = 1:numWins
                        %off set included
                        startWinPt = round(1+(winLen*(n-1)*fs));
                        endWinPt = round(min(winLen*n*fs,size(blockData,1)));
                        %Since all Nan set to 0 look for 0's instead of NaN
                        tmpNan(n,:) = sum(isnan(blockNan(startWinPt:endWinPt,:)),1);
                        
                        idxN = find(isnan(blockNan(startWinPt:endWinPt,1)) == 0);  %find all indices without Nan
                        x = [0; cumsum(diff(idxN)~=1)];                 %find consecutive non Nan values
                        
                        %check to see if there are any outages if there are
                        %not then proceed normally else do computation to
                        %componsate for Nans
                        if max(x) == 0
                            y = blockData(startWinPt:endWinPt,1:nChan);
                            R = corrcoef(y);
                            
                            %temporary soln to removing Nans that result
                            %from signal having 0 variance... why is that
                            %even happening???
                            R(isnan(R)) = 0;
                            
                            tmpEig(n,:) = eig(R)';
                            
                        else
                            %number of consecutive outs
                       	    numOut = max(x)+1;                                  
                            idxStart = ones(numOut,1);    %start ind                
                            idxEnd = ones(numOut,1);      %end ind
                            idxEnd(end) = length(x);      %last end point is end

                            %find all start and end index of outage
                            idxEnd(1:end-1) = find(diff(x) > 0);
                            idxStart(2:end) = find(diff(x) > 0) + 1;

                            %Start/end/and length of outage
                            startT = idxN(idxStart);   %values are not really time
                            endT = idxN(idxEnd);
                            outLen = (endT - startT);
                        
                            k = 1;
                            winEig = [];
                            for goodWin = 1:numOut
                                if outLen(goodWin) > 5                
                                
                                tmp = blockData(startT(goodWin):endT(goodWin),1:nChan);
                                R = corrcoef(tmp);
                                R(isnan(R)) = 0;
                                winEig(k,:) = eig(R)';
                                k = k + 1;
                                end
                            end
                            tmpEig(n,:) = mean(winEig,1);
                            
                        end
                end
                
            end
            
            eigVals{j} = tmpEig;
            numNan{j} = tmpNan;
            percentDone = 100 * j / numBlocks;
            msg = sprintf('Percent done: %3.1f',percentDone); %Don't forget this semicolon
            fprintf([reverseStr, msg])
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    fprintf('\n');
    eigVals = cell2mat(eigVals);
    numNan = cell2mat(numNan);
    save([datasetFN '_' outLabel '.mat'],'eigVals','-v7.3');
    end
end

end
