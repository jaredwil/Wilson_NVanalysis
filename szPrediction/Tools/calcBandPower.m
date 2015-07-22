function [alpha,beta,gamma,theta,delta,numNan] = calcBandPower(datasets,channels,band,winLen,outLabel,runIndex,blockLenSecs,parFlag,varargin)
%  [feat, numNan] = calcFeature_wil(datasets,channels,feature,winLen,outLabel,runIndex,blockLenSecs,parFlag,varargin)
%       This function is used for feature extraction over long data sets
%       able to select parflag to run processes in paralell.
%
% Author: Jared Wilson
% 4/8/15 v1 working
%
% NOTE: This function is a derivative of calcFeature function developed by
% Hoameng Ung
% Usage: calcFeature(datasets,channels,feature,winLen,outLabel,filtFlag,varargin)
% This function will divide IEEGDataset channels into blocks of
% and within these blocks further divide into winLen. Features
% will be calculated for each winLen and saved in a .mat matrix.
% Features calculated: power, LL, DCN
%
% Parallel Processing Usage: if parFlag set multiple sessions are
% established to pull data blocks simultaneosuly
%
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
%%%%%%%%%%%%%%%%%%%%%%%NOT REALLY USED%%%%%%%%%%%%%%%%%%%%
%   'band'         :  binary String input to define band of interest. Supported
%   bands in order are:  abgtd (xxxxx)  
%               alpha: 8-12 Hz
%                beta: 12-27 Hz
%               gamma: 27-45 Hz
%               theta: 3-8 Hz
%               delta: 0.2-3 Hz
%
%Output:
%   'feat'          :   Calculated feature for each window
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
CalcNumWins = @(xLen, fs, winLen, winDisp)floor((xLen-(winLen-winDisp)*fs)/(winDisp*fs));
DCNCalc = @(data) (1+(cond(data)-1)/size(data,2)); % DCN feature
AreaFn = @(x) mean(abs(x),1);
EnergyFn = @(x) mean(x.^2,1);
ZCFn = @(x) sum((x(1:end-1,:)>repmat(mean(x),size(x,1)-1,1)) & x(2:end,:)<repmat(mean(x),size(x,1)-1,1) | (x(1:end-1,:)<repmat(mean(x),size(x,1)-1,1) & x(2:end,:)>repmat(mean(x),size(x,1)-1,1)));
LLFn = @(x) mean(abs(diff(x)));
LLFn2 = @(X, winLen) conv2(abs(diff(X,1)),  repmat(1/winLen,winLen,1),'same');

%Important Const.
%freq Ranges of interest
alphaB = [8 12];
betaB = [12 27];
gammaB = [27 45];
thetaB = [3 8];
% delta = [0.2 3];
deltaB = [0 4];


%% Initialization
band = lower(band);
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
        %%%%%ESTABLISH MULTIPLE SESSIONS HERE!!!!!!!!!!!!%%%%%%%%%%%%%
        numPar = 6;     %This value is currenlty defined should be modifiable
        username = 'jaredwil'; % <<< These values should also be inputs somehow???
        pswd = 'jar_ieeglogin.bin';% <<<<<<<

        %determine number of blocks per parallel session
        parBlocks = floor(numBlocks/numPar);
        
        %off set index for each processor
        parOffset = parBlocks*blockLenSecs*fs;
        parFeat = cell(numPar,1);
        
        parAlpha = cell(numPar,1);
        parBeta = cell(numPar,1);
        parGamma = cell(numPar,1);
        parTheta = cell(numPar,1);
        parDelta = cell(numPar,1);
        
        parNan = cell(numPar,1);
        
        %%%REC CHANGE
        featTot = cell(parBlocks*numPar,1);
        
        featAlpha = cell(parBlocks*numPar,1);
        featBeta = cell(parBlocks*numPar,1);
        featGamma = cell(parBlocks*numPar,1);
        featTheta = cell(parBlocks*numPar,1);
        featDelta = cell(parBlocks*numPar,1);

        numNanTot = cell(parBlocks*numPar,1);
        
        parfor p = 1:numPar
                %warning('off')
                %Initilize variable to be used specifically in each proc.
                reverseStr = '';
                parfeat = cell(parBlocks,1);
                
                paralpha = cell(numPar,1);
                parbeta = cell(numPar,1);
                pargamma = cell(numPar,1);
                partheta = cell(numPar,1);
                pardelta = cell(numPar,1);
                
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
                    tmpFeat = zeros(numWins,nChan);
                    
                    tmpAlpha = zeros(numWins,nChan);
                    tmpBeta = zeros(numWins,nChan);
                    tmpGamma = zeros(numWins,nChan);
                    tmpTheta = zeros(numWins,nChan);
                    tmpDelta = zeros(numWins,nChan);
                    
                    tmpNan  = ones(numWins,nChan)*(winLen*fs);
                    
                else       %do normal feature extraction
                    blockNan = blockData; %copy block data with Nans 
                    blockData(isnan(blockData)) = 0; %NaN's turned to zero
                    nChan = numel(channels);
                    %calculate feature every winLen secs
                    numWins = ceil(size(blockData,1)/(winLen*fs));
                    tmpFeat = zeros(numWins,nChan);
                    
                    tmpAlpha = zeros(numWins,nChan);
                    tmpBeta = zeros(numWins,nChan);
                    tmpGamma = zeros(numWins,nChan);
                    tmpTheta = zeros(numWins,nChan);
                    tmpDelta = zeros(numWins,nChan);
                    
                    tmpNan  = zeros(numWins,nChan);
                    for n = 1:numWins
                        %off set included
                        startWinPt = round(1+(winLen*(n-1)*fs));
                        endWinPt = round(min(winLen*n*fs,size(blockData,1)));
                        %Since all Nan set to 0 look for 0's instead of NaN
                        tmpNan(n,:) = sum(isnan(blockNan(startWinPt:endWinPt,:)),1);
                        totSamp = winLen*fs;
                        
                        for c = 1:nChan
                        y = blockData(startWinPt:endWinPt,c);            
                        [PSD,F]  = pwelch(y,ones(length(y),1),0,length(y),fs,'psd');
                            
%                                 if(band(1) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpAlpha(n,c) = bandpower(PSD,F,alphaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpAlpha(n,c) = bandpower(PSD,F,alphaB,'psd');
                                        end                                     
%                                 elseif(band(2) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpBeta(n,c) = bandpower(PSD,F,betaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpBeta(n,c) = bandpower(PSD,F,betaB,'psd');
                                        end   
%                                 elseif(band(3) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpGamma(n,c) = bandpower(PSD,F,gammaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpGamma(n,c) = bandpower(PSD,F,gammaB,'psd');
                                        end 
%                                 elseif(band(4) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpTheta(n,c) = bandpower(PSD,F,thetaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpTheta(n,c) = bandpower(PSD,F,thetaB,'psd');
                                        end 
%                                 elseif(band(5) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpDelta(n,c) = bandpower(PSD,F,deltaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpDelta(n,c) = bandpower(PSD,F,deltaB,'psd');
                                        end 
%                                 end
                        end
                    end

                end

                parfeat{j} = tmpFeat;
                
                paralpha{j} = tmpAlpha;
                parbeta{j}  = tmpBeta;
                pargamma{j} = tmpGamma;
                partheta{j} = tmpTheta;
                pardelta{j} = tmpDelta;
                  
                parnumNan{j} = tmpNan;
            end
            
        parFeat{p} = parfeat; %Store the result from isolated process in 
        
        parAlpha{p} = paralpha;
        parBeta{p}  = parbeta;
        parGamma{p} = pargamma;
        parTheta{p} = partheta;
        parDelta{p} = pardelta;
        
        parNan{p} = parnumNan; %its appropriate location in cell
        end
    %put all cells back together in matrix
      feat = [];
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
       
       featAlpha(parBlocks*(tm-1)+1 : parBlocks*tm) = parAlpha{tm};
       featBeta(parBlocks*(tm-1)+1 : parBlocks*tm)  = parBeta{tm};
       featGamma(parBlocks*(tm-1)+1 : parBlocks*tm) = parGamma{tm};
       featTheta(parBlocks*(tm-1)+1 : parBlocks*tm) = parTheta{tm};
       featDelta(parBlocks*(tm-1)+1 : parBlocks*tm) = parDelta{tm};
       
       numNanTot(parBlocks*(tm-1)+1 : parBlocks*tm) = parNan{tm}; 
       
    end
    
    fprintf('\n');
    feat = cell2mat(featTot);
    
    alpha = cell2mat(featAlpha);
    beta  = cell2mat(featBeta);
    gamma = cell2mat(featGamma);
    theta = cell2mat(featTheta);
    delta = cell2mat(featDelta);
    
    numNan = cell2mat(numNanTot);

    save([datasetFN '_' outLabel '.mat'],'feat','-v7.3');
    
    save([datasetFN '_alpha_' outLabel '.mat'],'alpha','-v7.3');
    save([datasetFN '_beta_' outLabel '.mat'],'beta','-v7.3');
    save([datasetFN '_gamma_' outLabel '.mat'],'gamma','-v7.3');
    save([datasetFN '_theta_' outLabel '.mat'],'theta','-v7.3');
    save([datasetFN '_delta_' outLabel '.mat'],'delta','-v7.3');
    
    
    
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
        feat = cell(numBlocks,1);
        
        alpha = cell(numBlocks,1);
        beta = cell(numBlocks,1);
        gamma = cell(numBlocks,1);
        theta = cell(numBlocks,1);
        delta = cell(numBlocks,1);
        
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
                tmpFeat = zeros(numWins,nChan);
                
                tmpAlpha = zeros(numWins,nChan);
                tmpBeta = zeros(numWins,nChan);
                tmpGamma = zeros(numWins,nChan);
                tmpTheta = zeros(numWins,nChan);
                tmpDelta = zeros(numWins,nChan);
                
                
                tmpNan  = ones(numWins,nChan)*(winLen*fs);
            else       %do normal feature extraction
                blockNan = blockData; %copy block data with Nans 
                blockData(isnan(blockData)) = 0; %NaN's turned to zero
                nChan = numel(channels);
                %calculate feature every winLen secs
                numWins = ceil(size(blockData,1)/(winLen*fs));
                tmpFeat = zeros(numWins,nChan);
                
                tmpAlpha = zeros(numWins,nChan);
                tmpBeta = zeros(numWins,nChan);
                tmpGamma = zeros(numWins,nChan);
                tmpTheta = zeros(numWins,nChan);
                tmpDelta = zeros(numWins,nChan);
                
                tmpNan  = zeros(numWins,nChan);
                for n = 1:numWins
                    %off set included
                    startWinPt = round(1+(winLen*(n-1)*fs));
                    endWinPt = round(min(winLen*n*fs,size(blockData,1)));
                    %Since all Nan set to 0 look for 0's instead of NaN
                    tmpNan(n,:) = sum(isnan(blockNan(startWinPt:endWinPt,:)),1);
                    
                    totSamp = winLen*fs;
                    for c = 1:nChan
                        y = blockData(startWinPt:endWinPt,c);            
                        [PSD,F]  = pwelch(y,ones(length(y),1),0,length(y),fs,'psd');
                            
%                                 if(band(1) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpAlpha(n,c) = bandpower(PSD,F,alphaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpAlpha(n,c) = bandpower(PSD,F,alphaB,'psd');
                                        end                                     
%                                 elseif(band(2) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpBeta(n,c) = bandpower(PSD,F,betaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpBeta(n,c) = bandpower(PSD,F,betaB,'psd');
                                        end   
%                                 elseif(band(3) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpGamma(n,c) = bandpower(PSD,F,gammaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpGamma(n,c) = bandpower(PSD,F,gammaB,'psd');
                                        end 
%                                 elseif(band(4) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpTheta(n,c) = bandpower(PSD,F,thetaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpTheta(n,c) = bandpower(PSD,F,thetaB,'psd');
                                        end 
%                                 elseif(band(5) == '1')
                                        %scale
                                        if((1 - tmpNan(n,c)/totSamp) ~= 0)
                                        tmpDelta(n,c) = bandpower(PSD,F,deltaB,'psd')/(1 - tmpNan(n,c)/totSamp);
                                        else
                                        tmpDelta(n,c) = bandpower(PSD,F,deltaB,'psd');
                                        end 
%                                 end
                    end
                end
                
            end
            
            feat{j} = tmpFeat;
            
            alpha{j} = tmpAlpha;
            beta{j}  = tmpBeta;
            gamma{j} = tmpGamma;
            theta{j} = tmpTheta;
            delta{j} = tmpDelta;
            
            numNan{j} = tmpNan;
            percentDone = 100 * j / numBlocks;
            msg = sprintf('Percent done: %3.1f',percentDone); %Don't forget this semicolon
            fprintf([reverseStr, msg])
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    fprintf('\n');
    feat = cell2mat(feat);
    
    alpha = cell2mat(alpha);
    beta  = cell2mat(beta);
    gamma = cell2mat(gamma);
    theta = cell2mat(theta);
    delta = cell2mat(delta);
    
    numNan = cell2mat(numNan);

    %save([datasetFN '_' outLabel '.mat'],'feat','-v7.3');
    
    save([datasetFN '_alpha_' outLabel '.mat'],'alpha','-v7.3');
    save([datasetFN '_beta_' outLabel '.mat'], 'beta','-v7.3');
    save([datasetFN '_gamma_' outLabel '.mat'],'gamma','-v7.3');
    save([datasetFN '_theta_' outLabel '.mat'],'theta','-v7.3');
    save([datasetFN '_delta_' outLabel '.mat'],'delta','-v7.3');
    
    end
end

end
