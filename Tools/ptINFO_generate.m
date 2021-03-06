%Get statistics about pts. and generate excel
%Record Length -- rounded(days)
%number of CCS seizures per paitient

%Jared Wilson
%8/28/2015
% This is now a script that finds the sz annotations and return a matrix to
% be reference with the class of sz. Below is the legend for the type.
%    
%   0 - UCS (other electrical abnormality)
%   1 - CCS (clinically confirmed)
%   2 - CES (electricaly confirmed Sz)????
%   3 - NCS (not reported clinical Sz)
%
%%
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp

%%
% Define algorithm specifics 
usernm = 'jaredwil';
pswdBin = 'jar_ieeglogin.bin';
trPct = 0.7;
winLen = 30;
winDisp = 30;
szHorizon = 1; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

%attributes for ccs
numCCS = cell(1,numel(pt));
lengthRec = cell(1,numel(pt));
heading = {'Paitient' ;'Length of Recording (days)' ;'Number of CCS'};
%%
%begin function
%loop through all pts
for i = 1:numel(pt);  %%%%%TEMPORARY for debug

session = IEEGSession(pt{i},usernm,pswdBin) ;
fs = session.data.sampleRate;               %Find sampling Rate

%check to see if current pt. has annotations
if (size(session.data.annLayer,1) ~= 0)
    for layer = 1:size(session.data.annLayer,2)
        layerName = session.data.annLayer(layer).name;
        if(isempty(strfind(lower(layerName),'ccs')))
            continue;
        end

        %need to add to do this for only sz annot layers
        [~, timeUSec, ~] = getAllAnnots(session.data,layerName);

    end

    numCCS{i} = length(timeUSec);
    lengthRec{i} = round((session.data.channels(1).get_tsdetails.getDuration*(10^-6))/60/60/24);
    
else
    %no annotations
    disp(['WARNING: NO Annotations for this Pt. ' ... 
        'startT/endT returned as an empty matrix'])
    numCCS{i} = 0;
    lengthRec{i} = round((session.data.channels(1).get_tsdetails.getDuration*(10^-6))/60/60/24);

end

end


X = [pt; lengthRec; numCCS];
X = [heading X]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now write this CSV
fileName = 'C:\Users\Jared\Dropbox\Thesis stuff\ptInfo.csv';
  fid = fopen(fileName,'w');
%Print Heading  
  fprintf(fid,'%s, %s, %s\n',X{1,:});
  %Print each line independently
  for p = 1:numel(pt)
    fprintf(fid,'%s, %f, %f\n',X{p+1,:});
  end
  fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
