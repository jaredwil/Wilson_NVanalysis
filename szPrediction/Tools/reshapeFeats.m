%RESHAPE FEATURES %Could be a function????

% %%%%%%%%%%%%%%%%%%%%%%%RESHAPE FEATURES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1:2; %number of history samples
% samples = size(allFeat,1);
% features = size(allFeat,2);
% startSAMP = max(N)+1;
% M = (samples-length(N));  %number of time bins
% %create 'R' matrix for linear regression algorithm
% r = zeros(M, features*length(N)+1);
% for resIdx = 1:M
%     temp = allFeat(startSAMP + (resIdx-1) - N,:);   %temp is a temporary matrix    
%     r(resIdx,:) = [1 temp(:)'];
% end
% 
% allLabs = allLabs(startSAMP:end,:);
% allSzLab = allSzLab(startSAMP:end,:);
% 
% [~,IdxRemove] = findpeaks(allLabs);
% 
% r(IdxRemove,:) = [];
% allSzLab(IdxRemove,:) = [];
% allLabs(IdxRemove,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%