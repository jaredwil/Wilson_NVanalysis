%Long Term Analysis Test 1 All PT
clear all;
close all;
clc;

%Load .mat containing features
load('test1_info.mat')
energy = info{1};
ll = info{2};
numNan = info{3};
devEn = info{4};
devLL = info{5};

%find total number of NaNs in the first 60 days of each Pt.
for i = 1:length(numNan)
   tmp = numNan{i};   
   totNan(i) = sum(tmp);
end

%plot Win The number of averaged bins to be created for the plot 
%    60 - each bin is equal to 1 day
%    120 - each bin is equal to half a day
plotWin = 60;
%%%%%%%%%%%%%%%%%%%%%%%%%%Line Length Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%
dayLL = zeros(size(ll,1),plotWin);

%Loop through all Pt.
for i = 1:length(ll)
    %retrieve ll from current pt. 
    tmpLL = ll{i};
    
    %reshape into matrix where each row represents 1 day
    tmpLL = reshape(tmpLL, length(tmpLL)/plotWin, plotWin);
    
    %loop through all days
    for j = 1:size(tmpLL,2)
        
       %remove all zero values from that day then average this prevents the
       %zero values from scewing the average
       tmp = tmpLL(:,j);
       tmp(tmp == 0) = [];
       
       if tmp %only calculate the average if tmp is not empty
           dayLL(i,j) = mean(tmp);
           cellLL{j} = tmp;
       end
        
    end
end

%remove outliers
skipInd = [1:9 12:14];
dayLL = dayLL(skipInd,:);

%average each day across all pt.
for i = 1:size(dayLL,2)
    %All zero feature values removed from that day
    tmpDay = dayLL(:,i);
    tmpDay(tmpDay == 0) = [];
    
    %find std and mean using all pt with non zero features
    llAv(i) = mean(tmpDay);
    llStd(i) = std(tmpDay);
end

% llAv = llAv(5:end);
% llStd = llStd(5:end);


figure(1)
% plot(dayEnergy')
errorbar(llAv, llStd)
xlabel('Day')
ylabel('Line Length')
title('Average Time Windowed Line Length Across All Patients')

figure(2)
plot(llAv)

figure(3)
plot(dayLL')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Energy Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Energy analysis done in the same way as above

dayEnergy = zeros(size(energy,1),plotWin);
%Pt 5 is a good one.
for i = 1:length(energy)
    tmpEnergy = energy{i};
    
    tmpEnergyNa = tmpEnergy;
    tmpEnergyNa(tmpEnergy == 0) = [];
 
    tmpEnergy = reshape(tmpEnergy, length(tmpEnergy)/plotWin, plotWin);
    
    for j = 1:size(tmpEnergy,2)
       tmp = tmpEnergy(:,j);
       tmp(tmp == 0) = [];
       
       if tmp
           dayEnergy(i,j) = mean(tmp);
           cellEng{j} = tmp;
       end
        
    end

end

for i = 1:size(dayEnergy,2)
    tmpDay = dayEnergy(:,i);
    
    tmpDay(tmpDay == 0) = [];
        
    engAv(i) = mean(tmpDay);
    engStd(i) = std(tmpDay);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%IMPORTANT GRAPH%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(33)
% plot(dayEnergy')
errorbar(engAv*10^-3, engStd*10^-3)
xlabel('Day')
ylabel('Energy [nJ]')
title('Average Time Windowed Energy Across All Patients')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(34)
plot(engAv)

figure(35)
plot(dayEnergy')