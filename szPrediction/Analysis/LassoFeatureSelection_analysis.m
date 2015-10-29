%Lasso Feature vector explained
%this script will disect the feature coef returned by lasso to see which
%features have been zeroed
numCh = 16;



%%
% Define algorithm specifics 
usernm = 'jaredwil';
pswdBin = 'jar_ieeglogin.bin';
trPct = 0.7;
winLen = 30;
winDisp = 30;
szHorizon = 2; %hours

% patients of interest on ieeg portal
pt = {'NVC1001_25_001' 'NVC1001_25_002' 'NVC1001_25_004' ...
    'NVC1001_25_005' 'NVC1001_24_001' 'NVC1001_24_002' 'NVC1001_24_004' ...
    'NVC1001_24_005' 'NVC1001_23_002' 'NVC1001_23_003' 'NVC1001_23_004' ...
    'NVC1001_23_005' 'NVC1001_23_006' 'NVC1001_23_007'};

featInfo_all = zeros(numel(pt),14);
featInfo_ind = zeros(numel(pt),14);
chInfo_all   = zeros(numel(pt),16);

sumBetas = zeros(numel(pt),14);
percBetas = zeros(numel(pt),14);

%%
%begin function
%loop through all pts
for i = 1:numel(pt);  

    
%     featLab = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Res LassoLib\' pt{i} '_bestLasso5.mat'];
    featLab = ['H:\jaredwil\Lasso Results\Res 8-29\' pt{i} '_bestLasso5.mat'];

    try   
        load(featLab);
        
    catch
        disp('This pt. does not have a save .mat to load from')
        continue;
    end
    
    coef = lassoRes.coef.totmin;
    int = lassoRes.int(3);

    %CURRENT FEATURE VECTOR FORMAT AS OF 8/17/2015 check to make sure curernt
    % feats = [Amp LL nlEng alpha beta delta gamma theta fxskew fxkurt psdCorr fycorr fycov HFD];

    amp    = coef(1:16);
    ll     = coef(1*16+1:16*2);
    nlEng  = coef(2*16+1:16*3);
    alpha  = coef(3*16+1:16*4);
    beta   = coef(4*16+1:16*5);
    delta  = coef(5*16+1:16*6);
    gamma  = coef(6*16+1:16*7);
    theta  = coef(7*16+1:16*8);
    fxskew = coef(8*16+1:16*9);
    fxkurt = coef(9*16+1:16*10);
    psdCorr = coef(10*16+1:10*16+120);
    fycorr  = coef(10*16+120+1:10*16+240);
    fycov   = coef(10*16+240+1:10*16+360);
    HFD     = coef(10*16+360+1:10*16+360+16);

    %find how many channels were zerod for each feature
    numNotZero = [nnz(amp) nnz(ll) nnz(nlEng) nnz(alpha) nnz(beta) nnz(delta) ...
        nnz(gamma) nnz(theta) nnz(fxskew) nnz(fxkurt) nnz(psdCorr) nnz(fycorr) ...
        nnz(fycov) nnz(HFD)];
    
    totNZ(i) = sum(numNotZero);
    
    sumBetas(i,:) = [sum(abs(amp)) sum(abs(ll)) sum(abs(nlEng)) sum(abs(alpha)) sum(abs(beta)) sum(abs(delta)) ...
        sum(abs(gamma)) sum(abs(theta)) sum(abs(fxskew)) sum(abs(fxkurt)) sum(abs(psdCorr)) sum(abs(fycorr)) ...
        sum(abs(fycov)) sum(abs(HFD))];
    
    percBetas(i,:) = abs(sumBetas(i,:))./sum(abs(sumBetas(i,:)));
    
    featLens = [16 16 16 16 16 16 16 16 16 16 120 120 120 16];
    featInfo_ind(i,:) = bsxfun(@rdivide, numNotZero, featLens);
    featInfo_all(i,:) = numNotZero./sum(numNotZero);


    %is there a way to find out which channels contain important information
    %for regression????
    %For ch-by-ch features this is very easy... for the cov and corr matrices
    %this is a bit more challenging.
    rowCh = [];
    colCh = [];
    for j = 1:15
        rowCh = [rowCh 1:j];
        colCh = [colCh ones(1,j)*(j+1)];
    end

    totCh = rowCh | colCh;

    %go through all the channels
    chNotZero = zeros(1,16);
    for ch = 1:numCh
        chIdx = (rowCh == ch | colCh == ch);

        chFeat = [amp(ch) ll(ch) nlEng(ch) alpha(ch) beta(ch) delta(ch) ...
            gamma(ch) theta(ch) fxskew(ch) fxkurt(ch) psdCorr(chIdx)' fycorr(chIdx)' ...
            fycov(chIdx)' HFD(ch)];
        chNotZero(ch) = nnz(chFeat);
    end
    
    chInfo_all(i,:) = chNotZero./sum(numNotZero);


end

figure(1)
plotfeatInfo(sumBetas);
title('Sum of Weights Associated with Each Feature')
% set(gcf,'Position',get(0,'Screensize')); 
% plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\Meeting Results 9-20\featInfo'];
plotName = ['H:\jaredwil\Lasso Results\featureAnalysis\sumWeights'];
saveas(gcf,plotName,'jpg')
set(gca,'FontSize',10);
set(gca,'LineWidth',2);

figure(2)
plotfeatInfo(percBetas);
title('Precent of Total Weights Associated with Each Feature')
colorbar();
% set(gcf,'Position',get(0,'Screensize')); 
set(gca,'FontSize',13);
set(gca,'LineWidth',2);
% plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\Meeting Results 9-20\featInfo2'];
plotName = ['H:\jaredwil\Lasso Results\featureAnalysis\percWeights'];
saveas(gcf,plotName,'jpg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%OLD PLOTS
%These plots were very confusing



%plot some stuff
% figure(1)
% plotfeatInfo(featInfo_all*100);
% title('Precent of Weights Associated with Feature')
% % set(gcf,'Position',get(0,'Screensize')); 
% plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\Meeting Results 9-20\featInfo'];
% saveas(gcf,plotName,'jpg')
% set(gca,'FontSize',10);
% set(gca,'LineWidth',2);
% 
% figure(2)
% plotfeatInfo(featInfo_ind*100);
% title('Precent of Non-Zero Weight for every Feature ')
% % set(gcf,'Position',get(0,'Screensize')); 
% % set(gca,'FontSize',15);
% % set(gca,'LineWidth',2);
% plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\Meeting Results 9-20\featInfo2'];
% saveas(gcf,plotName,'jpg')
% set(gca,'FontSize',10);
% set(gca,'LineWidth',2);
% 
% figure(3)
% plotchannelInfo(chInfo_all*100);
% title('Percent of Weights Associated with Each Channel')
% % set(gcf,'Position',get(0,'Screensize')); 
% % set(gca,'FontSize',15);
% % set(gca,'LineWidth',2);
% plotName = ['C:\Users\Jared\Dropbox\NVanalysis_data\Meeting Results 9-20\chInfo'];
% saveas(gcf,plotName,'jpg')
% set(gca,'FontSize',10);
% set(gca,'LineWidth',2);