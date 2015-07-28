function [ trainFeats, trainLabels] = szPred_train(pt, usernm, pswdBin, trPct, winLen, winDisp, szHorizon  )
%[ trainFeats, trainLabels ] = szPred_train(pt, usernm, pswdBin, trPct, winLen, winDisp, szHorizon  )
%   szPred_train takes the current pt and retrieves the trainFeats and
%   trainLabels based on the inputs provided by the user. The szHorizon
%   defines the 'preictal' phase where features will be extracted from all
%   lables that are larger than the defined szHorizon are interictal
%   feature vectors. 
%
%Inputs:
%       pt - string containg pt. of interest on ieeg portal
%       usernm - user name for ieeg
%       pswdBin - .bin file for ieeg accout (make sure is in path
%       somewhere)
%       trPct   - the percent of data to be trained on
%       winLen  - window length for feature extraction
%       winDisp - displacement of consecutive windows
%       szHorizon - defines what labels will determine as 'preictal' (refer
%       to Liturature for suggested values range: 30min - 2 hrs). 
%
%Output:
%      trainFeats - M x N matrix where M is the number of observations and
%      N is the number of features extracted from each time window defined
%      by winLen and winDisp. (NOT NORMALIZED)
%
%      trainLabels - M x 1 matrix containing the corrosponding lables to
%      observations in trainFeats. Labels are time to sz in minuites.
%
%%

session = IEEGSession(pt,usernm,pswdBin) ;
fs = session.data.sampleRate;               %Find sampling Rate
 
%get ALL seziure start and end times from anotations
[startT, endT] = getSzTime(session);

numTr = floor(length(startT)*trPct);  %number of training sz

%split the sz into train/test sets
trainST = startT(1:numTr);
trainET = endT(1:numTr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%TEMPORARY TO INCREASE SPEED!!!!!!
%if more than 30 training sz then only use first 30 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get training Data
%Sz horizon data
[train, predIdx] = getSzFeats(session, trainST, trainET, szHorizon, winLen, winDisp);
%interictal data retrieved from a random collection of samples from various
%hours in the training data.

Nlim = size(train,1);
[trainInt] = InterictRand_train( session, predIdx, winLen, winDisp, Nlim);

%isolate time lables and training feats
szhorzTrain = train(:,2:end);
szhorzLabels = train(:,1);
intLables = ones(size(trainInt,1),1)*(szHorizon*60*60 + winLen);

trainLabels = [szhorzLabels; intLables];
trainFeats = [szhorzTrain;trainInt];

end

