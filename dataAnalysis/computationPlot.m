close all;
clc;

% addpath(genpath('NVanalysis_data'))

load('timeTest.mat')
% par12 = load('par12.mat');
% par12 = par12.timePar;
par9 = load('par9.mat');
par9 = par9.timePar;
par3 = load('par3.mat');
par3 = par3.timePar;

tSer = f{1};
tPar = f{2};
check = f{3};
testN = [1 2 3 4 5 10 15 20 25 30 35 40 45 50 60 70 80 90 100 125 150 200];

tDif = tSer - tPar;

figure(1)
% set(gcf,'Position',get(0,'Screensize')); 
%make background white
set(gcf,'Color','w');
plot(testN,tSer,'k.-','MarkerSize',15,'LineWidth',1.5)
hold all;
grid on;
plot(testN,par3,'g.-','MarkerSize',15,'LineWidth',1.5)
plot(testN,tPar,'b*-','MarkerSize',8,'LineWidth',1.5)
plot(testN,par9,'m.-','MarkerSize',15,'LineWidth',1.5)
% plot(testN,tDif,'r.--','MarkerSize',10)
% plot(testN,par12,'k.-','MarkerSize',10)
% plot(testN,par9,'m.-','MarkerSize',10)

xlabel('Number of Days')
ylabel('Computation Time (sec)')
legend('Serial','Parallel (3 CPUs)','Parallel (6 CPUs)','Parallel (9 CPUs)','location','best')
set(gca,'FontSize',18);
set(gca,'LineWidth',2);