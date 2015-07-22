function [ startT, endT ] = getSzTime(ptSession)
%[startT] = getSzStartT(ptSession) 
%   Get ALL seizure start times from the pt of interest. The seizure start
%   times are extracted from the annotation layers. All Layers are
%   extracted.
%
%Inputs:
%       ptSession - initiated iEEG session
%
%Output:
%      startT - start time of all seizures in seconds
%%

timeUSec = cell(size(ptSession.data.annLayer,2),1);
startT = [];
endT = [];
%check to see if current pt. has annotations
if (size(ptSession.data.annLayer,1) ~= 0)
    for layer = 1:size(ptSession.data.annLayer,2)
        layerName = ptSession.data.annLayer(layer).name;

        %need to add to do this for only sz annot layers
        [~, timeUSec{layer}, ~] = getAllAnnots(ptSession.data,layerName);

        startT = [startT; (timeUSec{layer}(:,1)*1e-6)];
        endT = [endT; (timeUSec{layer}(:,2)*1e-6)];
 
    end

else
    %no annotations
    disp(['WARNING: NO Annotations for this Pt. ' ... 
        'startT/endT returned as an empty matrix'])
    startT = [];
    endT = [];
end

%sort the start and end times from smalles to largest
startT = sort(startT);
endT = sort(endT);

end

