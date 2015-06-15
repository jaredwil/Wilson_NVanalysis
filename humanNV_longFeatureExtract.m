%This scrip is designed to run through the Human Neuro Vista recording and
%extract various time-windowed features for analysis. 



% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath('lib')
addpath('NV Analyze Toolbox')

%%
% <latex>
% \section*{1. Section 1}
% </latex>

%%
% <latex>
% \begin{enumerate}
%   \item Find sampling rate and spectral resolution of sample. 
% \end{enumerate}
% </latex>


session = IEEGSession('NVC1001_23_002','jaredwil','jar_ieeglogin.bin') ;
fs = session.data.sampleRate;               %Find sampling Rate


%get one seizure and view it
% feat1 = calcFeature_wil(session.data,1:16,'ll',2,'NVtest',[2931700 2931700+100], 10,  0);
% [feat2, numNan] = calcFeature_wil(session.data,1:16,'ll',0.5,'NVtest',[2931700 2931700+100], 10,  0);
% 
% plot(sum(feat1(:,1),2))
% plot(sum(feat2(:,1),2))

%%
%Get power and ll in one second windows set in 1 day increments
day = 86400; %sec
hour = 3600; %sec
min = 60; %sec;

% [ll, numNan] = calcFeature_wil(session.data,1:16,'ll',sec,'NVtest','all', hour,  0);
 [ll, numNan] = calcFeature_wil(session.data,1:16,'ll',30,'llTest2',[0 60*day], hour*2,  0);
 
 
ll = sum(ll,2)/size(ll,2);   %compute average ll over all 16 channels for each window
numNan = sum(numNan,2)/size(numNan,2); %computer average numNans over all 16 channes for each win

%%%%%%%%%%%%%%%%%%%%%%%s%%%%%%FIND FEATURE FOR EACH DAY CODE%%%%%%%%%%%%%%%

%find how many full days are in the data set;
% numDays = round(size(ll,1)/day);
numDays = round(size(ll,1)/(120));

%reshape the matrix so each column represents each day
llDays = reshape(ll(1:(numDays*120)),120,numDays);
%remove all 0 ll values which means that windows contained all Nans

%Compute the average and standard deviation in each day. Cell is used to
%contain each day because length varies dependent on how many Nan Windows
%exsit
llDays2 = {};
avgLL = zeros(numDays,1);
stdLL = zeros(numDays,1);
for i=1:size(llDays,2)
    llDays2{i} = llDays(:,i);
    llDays2{i}(llDays2{i} == 0) = [];
    avgLL(i) = mean(llDays2{i});
    stdLL(i) = std(llDays2{i});
end


%%%%%%%%%%%%%%%%%%%%%%HOURS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %find how many full days are in the data set;
% numhours = round(size(ll,1)/hour);
% %reshape the matrix so each column represents each day
% llhours = reshape(ll(1:(numhours)),hour,numhours);
% %remove all 0 ll values which means that windows contained all Nans
% 
% %Compute the average and standard deviation in each day. Cell is used to
% %contain each day because length varies dependent on how many Nan Windows
% %exsit
% llDays2 = {};
% avgLL = zeros(numDays,1);
% stdLL = zeros(numDays,1);
% for i=1:size(llDays,2)
%     llDays2{i} = llDays(:,i);
%     llDays2{i}(llDays2{i} == 0) = [];
%     avgLL(i) = mean(llDays2{i});
%     stdLL(i) = std(llDays2{i});
% end


h = plot(avgLL);

info = {};
info{1} = avgLL;
info{2} = stdLL;

saveas(h, 'avgLL.png','png')
save('info.mat','info')

% 
% power = calcFeature_wil(session.data,1:16,'power',86400,'NVtest','all', 86400,  0);
% 
% 
% figure(1)
% 
% subplot(211)
% plot(ll);
% xlabel('Day')
% ylabel('Average Line Length Over all Channels')
% subplot(212)
% plot(power);
% xlabel('Day')
% ylabel('Power over all Channels')









