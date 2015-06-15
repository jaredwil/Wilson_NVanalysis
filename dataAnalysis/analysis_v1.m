%Long Term Analysis Test 1 All PT
clear all;
close all;
clc;

load('test1_info.mat')


energy = info{1};
ll = info{2};
numNan = info{3};
devEn = info{4};
devLL = info{5};


for i = 1:length(numNan)
   tmp = numNan{i};   
   totNan(i) = sum(tmp);
end


plotWin = 60;

%%%%%%%%%%%%%%%%%%%%%%%%%%Line Length Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%
dayLL = zeros(size(ll,1),plotWin);
%Pt 5 is a good one.
for i = 1:length(ll)
    tmpLL = ll{i};
    
    tmpLLNa = tmpLL;
    tmpLLNa(tmpLL == 0) = [];
%     normAv = mean(tmpEnergyNa);
%     normStd = std(tmpEnergyNa);
%     
%     tmpEnergy = (tmpEnergy - normAv)./normStd;
    
    tmpLL = reshape(tmpLL, length(tmpLL)/plotWin, plotWin);
    
    for j = 1:size(tmpLL,2)
       tmp = tmpLL(:,j);
       tmp(tmp == 0) = [];
       
       if tmp
           dayLL(i,j) = mean(tmp);
           cellLL{j} = tmp;
       end
        
    end
%     dayEnergy(i,:) = mean(tmpEnergy,1);
%   plot(energy{i})
end

skipInd = [1:9 12:14];
dayLL = dayLL(skipInd,:);
for i = 1:size(dayLL,2)
    tmpDay = dayLL(:,i);
    
    tmpDay(tmpDay == 0) = [];
        
    llAv(i) = mean(tmpDay);
    llStd(i) = std(tmpDay);
      
%     engAv = mean(dayEnergy,1);
%     engStd = std(dayEnergy,1);

    
end

llAv = llAv(5:end);
llStd = llStd(5:end);


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
dayEnergy = zeros(size(energy,1),plotWin);
%Pt 5 is a good one.
for i = 1:length(energy)
    tmpEnergy = energy{i};
    
    tmpEnergyNa = tmpEnergy;
    tmpEnergyNa(tmpEnergy == 0) = [];
%     normAv = mean(tmpEnergyNa);
%     normStd = std(tmpEnergyNa);
%     
%     tmpEnergy = (tmpEnergy - normAv)./normStd;
    
    tmpEnergy = reshape(tmpEnergy, length(tmpEnergy)/plotWin, plotWin);
    
    for j = 1:size(tmpEnergy,2)
       tmp = tmpEnergy(:,j);
       tmp(tmp == 0) = [];
       
       if tmp
           dayEnergy(i,j) = mean(tmp);
           cellEng{j} = tmp;
       end
        
    end
%     dayEnergy(i,:) = mean(tmpEnergy,1);
%   plot(energy{i})
end

for i = 1:size(dayEnergy,2)
    tmpDay = dayEnergy(:,i);
    
    tmpDay(tmpDay == 0) = [];
        
    engAv(i) = mean(tmpDay);
    engStd(i) = std(tmpDay);
      
%     engAv = mean(dayEnergy,1);
%     engStd = std(dayEnergy,1);

    
end

engAv = engAv(5:end);
engStd = engStd(5:end);

% tmpEnergy = energy{5};
% tmpEnergy = reshape(tmpEnergy, 2880, 60);
% dayEnergy = mean(tmpEnergy,1);

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