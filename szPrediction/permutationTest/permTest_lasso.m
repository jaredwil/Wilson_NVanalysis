%Jared Wilson 
%Permutation TEST
%Creat permutation of labels and train a new model to see if original
%results are signif. (p < 0.05)

%%
% clear the workspace and console
clear all; close all; clc;
warning('off')
addpath(genpath('ieeg-matlab-1.8.3'))
addpath(genpath('Wilson_NVanalysis'))
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data'))  %this is where .mat file are contained on local comp
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\30secFeats'))  %this is where .mat file are contained on local comp
% addpath(genpath('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\5minFeats'))  %this is where .mat file are contained on local comp
addpath(genpath('H:\jaredwil\szPred_feats\Win_5min'))  %this is where .mat file are contained on local comp


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

%%%%%%%%%%%%%%%%%%%%%% start par pool%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Establish parpool
[p, poolsize ] = initParPool( 10 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lassoTime = zeros(1,numel(pt));
%%
%begin function
%loop through all pts
for i = 1:numel(pt);  %%%%%TEMPORARY for debug
     close all;
     clc;
%   label = [pt{i} '_szPred_30secFeats.mat'];
    label = [pt{i} '_szPred_5minFeats.mat'];

    tylabel = [pt{i} '_szTypeLabels5.mat'];

    try   
        load(label);
        load(tylabel);
    catch
        disp('This pt. does not have a save .mat to load from')
        continue;

    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% get feature data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    train   = data.train;
    test    = data.test;
    % avgFeats = data.mean;
    % stdFeats = data.std;

    %isolate labels and feats in training data
    trainLabels = train(:,1);
    trainFeats = train(:,2:end);
    testLabels = test(:,1);
    testFeats = test(:,2:end);

    %Remove intericatal Data from training data & Testing Data
    trainFeats(trainLabels > szHorizon*60*60,:) = [];
    trainLabels(trainLabels > szHorizon*60*60,:) = [];
    testFeats(testLabels > szHorizon*60*60,:) = [];
    testLabels(testLabels > szHorizon*60*60,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%% get sz time and type labels%%%%%%%%%%%%%%%%%%%%

    trainSzT = szType.labels.train;
    testSzT = szType.labels.test;
    
    szTimes = [szType.szT.train;szType.szT.test];
    
    szTimes(:,1) = szTimes(:,1)/60/60/24;
    
    avgFeats = mean(trainFeats,1);
    stdFeats = std(trainFeats,[],1);

    
    allFeat = [trainFeats; testFeats];
    allFeat = bsxfun(@rdivide, bsxfun(@minus,allFeat,avgFeats), stdFeats);
    
    trainFeats =  bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);
    testFeats =  bsxfun(@rdivide, bsxfun(@minus,testFeats,avgFeats), stdFeats);
    
    
    
    allLabs  = [trainLabels;testLabels];
    allSzLab  = [trainSzT;testSzT];
    
    
%%%%%%%%%%%%%% REMOVE SZ THAT ARE NOT CCS (CLINICALLY CONFIRMERD) %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ONLY LOOK AT CCS
    allFeat(allSzLab(:,2)~=1,:) = [];
    allLabs(allSzLab(:,2)~=1) = [];
    allSzLab(allSzLab(:,2)~=1,:) = [];

    %load feature weights
    % load('C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\resultsV1\NVC1001_23_005_bestLasso_CCS.mat')

    % featLab = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Res LassoLib\' pt{i} '_bestLasso5.mat'];
    featLab = ['H:\jaredwil\Lasso Results\Res 8-29\' pt{i} '_bestLasso5.mat'];

    load(featLab);


    
    fCorr = lassoRes.coef.corr;
    tIntCorr = lassoRes.int(1);
    tIntCorr2 = repmat(tIntCorr, size(allFeat,1),1);
    predCorr = allFeat*fCorr + tIntCorr2;

    
    
    fMin = lassoRes.coef.min;
    tIntMin = lassoRes.int(2);
    tIntMin2 = repmat(tIntMin, size(allFeat,1),1);
    predMin = allFeat*fMin + tIntMin2;


    f = lassoRes.coef.totmin;
    tInt = lassoRes.int(3);
    tInt2 = repmat(tInt, size(allFeat,1),1);
    testInt = repmat(tInt, size(testFeats,1),1);

    pred = allFeat*f + tInt2;
    predTest = testFeats*f + testInt;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Validate Model!
%this was to find appropriate lambda range
%     numLam = 75;
%     disp(['Training Lasso Model on Patient: ' pt{i} ' (NO INTERICTAL)'])
%     lambda   = linspace(1e4,1e6,numLam);
% %     lambda   = logspace(4,6,numLam);
%     numFolds = 10;
%     modCh = 'se';
%     trainFeats =  bsxfun(@rdivide, bsxfun(@minus,trainFeats,avgFeats), stdFeats);
% 
%     try
%         [ fre, tIntre, numFeats, solnLambda, lassoTime(i) ] = trainLassoLib( trainFeats,trainLabels,numFolds,lambda, modCh, pt{i});
%     catch
%         disp('Could Not Find Solution')
%     end
%     
%     tIntre2 = repmat(tIntre, size(allFeat,1),1);
%     predre = allFeat*fre + tIntre2;
% 
%     disp(['DONE Training Lasso Model on Patient: ' pt{i}])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%     %visual assesment
%     figure(1)
%     plot(allLabs);
%     hold on;
%     plot(smooth(pred,5));
%     plot(smooth(predCorr,5));
%     plot(smooth(predMin,5));
%     plot(smooth(predre,5));
%     legend('real','se','corr','min','re')




    nullN = 200;  %number of permutation that will be tested
    %Create a random label matrix mx1000 where m is the number of labels
    permLab_train = zeros(length(trainLabels),nullN);


    %average feature across sz do pca then cluster. 
    [numSz_train,~, uIdx] = unique(trainSzT(:,1),'stable'); 
    [numSz_test,~, uIdxTest] = unique(testSzT(:,1),'stable'); 
    
    timeDay   = zeros(numel(numSz_train),1);
    timeSec   = zeros(numel(numSz_train),1);
    winDur_sz = zeros(numel(numSz_train),1);
    for sz = 1:numel(numSz_train);
        
        winDur_sz(sz) = length(trainLabels(uIdx == sz));
        
        %create permutation labels rather than random labels
        for permIt = 1:nullN
            tmp = trainLabels(uIdx == sz);
            permord = randperm(winDur_sz(sz));
            permLab_train(uIdx == sz,permIt) = tmp(permord);
        end
        
        timeDay(sz) = mode(trainSzT(uIdx == sz,1))/60/60/24;
        timeSec(sz) = mode(trainSzT(uIdx == sz,1));
        
    end
        
    %eliminate sz that do not have many observations because
    %it messes up the null hypothesis
    % this does two things not really sure if I should be doing both at the
    % same time The first part eliminates short sz to remove false high sz
    % correlations, the second eliminates sz in the first 15 days. Features
    % Normalization
    idxElim = (winDur_sz < 6) | (timeDay < 15); 
    timeDay(idxElim)       = [];
%     timeSec(idxElim)       = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Permutation TEST 
    numLam = 5;
%     lambda   = linspace(1e4,1e6,numLam);
%     lambda   = linspace(5e3,1e5,numLam);

%     lambda   = [1000 5000 10000 50000 100000];
    lambda   = [8000 10000 12000 25000 100000];

%     permCorr = zeros(numel(numSz_train),numLam,nullN);  %correlation of permuted pred lables
%     permBeta = zeros(size(trainFeats,2),numLam,nullN); 
        
    permCorr = cell(nullN,1);  %correlation of permuted pred lables
    permBeta = cell(nullN,1);
    %Train Mode for each iteration
    parfor_progress(nullN);
    disp(['Training Pt.' pt{i}])
    
    tic;
    parfor permIt = 1:nullN
        
%         disp(['Creating Null Distribution Permutation ' num2str(permIt) '/' num2str(nullN)])
        
        w = zeros(size(trainFeats,2),numLam);
        for lam = 1:length(lambda)
            %change this functiion manually if you wish to use a dif method
            [w(:,lam), ~, ~] = LassoBlockCoordinate(trainFeats,permLab_train(:,permIt),lambda(lam));
%             parfor_progress;
%             disp(num2str(lam))
        end
        
        permInt = repmat(mean(permLab_train(:,permIt)),1,numLam) - mean((trainFeats*w),1);

        tInt2 = repmat(permInt, size(trainFeats,1),1);

        permPred = trainFeats*w + tInt2;

        tTest = repmat(testLabels,1,numLam);

%         tmp = corr(tTest,permPred);
%         tLasso_coor(cvIter,:) = tmp(1,:);
        itCorr = zeros(numel(numSz_test),numLam);  %correlation of permuted pred lables
        winDur_test = zeros(numel(numSz_test),1);
        testSz_type = zeros(numel(numSz_test),1);
        for sz = 1:numel(numSz_test);
            %permuted predicted labels for each sz rather than random
            %generation
%             permCorr(sz,:,permIt) = corr(tTest((uIdxTest == sz),:),permPred(uIdxTest == sz,:));
            
            winDur_test(sz) = length(testLabels(uIdxTest == sz));
            testSz_type(sz) = mode(testSzT(uIdxTest == sz,2));
            
            tmpCorr = corr(tTest((uIdxTest == sz),:),permPred(uIdxTest == sz,:));
            itCorr(sz,:) = tmpCorr(1,:);
        end
        idxElim = (winDur_test < 6) | (testSz_type ~= 1); %remove if short Sz or not CCS sz.

        %remove the short sz
        itCorr(idxElim,:)    = [];
        %save the beta values for the permutation
%         permBeta(:,:,permIt) = cell2mat(w);
        permCorr{permIt} = itCorr;
        permBeta{permIt} = w;
        
        
        parfor_progress;
    end
    permTime = toc;  %gives me an idea of how long thing take
    
    disp(['It took ' num2str(permTime) ' to compute the null distribution using ' ...
        num2str(nullN) ' different permutations.'])
    
    notElim_sz = size(permCorr{1},1);

    
    permCorr = cell2mat(permCorr);
    permBeta = cell2mat(permBeta);
    
    permCorr2 = reshape(permCorr,[notElim_sz, nullN,numLam]);
    permBeta2 = reshape(permBeta,[size(testFeats,2),nullN,numLam]);
    
%%%%%%compute test correlation to compare against null model


    testInt = repmat(tInt, size(testFeats,1),1);
    predTest = testFeats*f + testInt;
    testCorr = zeros(numel(numSz_test,1));
    winDur_test = zeros(numel(numSz_test),1);
    testSz_type = zeros(numel(numSz_test),1);
    for sz = 1:numel(numSz_test);
        winDur_test(sz) = length(testLabels(uIdxTest == sz));
        testSz_type(sz) = mode(testSzT(uIdxTest == sz,2));

        testCorr(sz,:) = corr(predTest(uIdxTest == sz),testLabels(uIdxTest == sz));
    end
    idxElim = (winDur_test < 6) | (testSz_type ~= 1); %remove if short Sz or not CCS sz.
    testCorr(idxElim) = [];
    
%%%%%%%%%%%%%%%%%%%create null distributions 
%%
%null distributions of values for each weight

%each cell contains a B x N matrix where N is the number of
%permutations/observations and B are the different weights

%each cell contains the beta distrubtion for each lambda value
% 
% fInterp = repmat(f,[1,nullN]);
% betaDist = cell(numLam,1);
% for d = 1:numLam
%    betaN = permBeta2(:,:,d) ;
%    betaDist{d} = betaN;
%    
%    %create plot
% %    betaRatio = fInterp./(betaN + 2e-9);  %add very small number to avoid dividing by 0
%    %plot histogram of the non-zero weights
%    figure(d+30)
%    [C1, C1x] = hist(betaN(betaN ~= 0),100);
%     C1 = C1./sum(C1);
%     bar(C1x,C1)
%     axis('tight')
%     title(['\lambda = ' num2str(lambda(d))])
%     hold on;
%     [C2, C2x] = hist((f~= 0),10);
%     C2 = C2./sum(C2);
%     bar(C2x,C2)
% end


%null distrubtion of correlation for each sz
%each cell contains a Sz x N matrix where N is the number of
%permutations/observations and Sz are individual Sz
lamdaLegend = strsplit(num2str(lambda));
szDist = cell(numLam,1);
probDen  = cell(notElim_sz,numLam);
corrNull = cell(notElim_sz,numLam);
yDist    = cell(notElim_sz,numLam);

for d = 1:numLam
    
   indszDist =  permCorr2(:,:,d); 
   szDist{d} = indszDist;
   

   %loop through all sz to create a plot of each... lots of plots I know
   for sz = 1:notElim_sz
       %fit the data using a kernel distribution Epanechnikov kernel
       %function
       probDen{sz,d} = fitdist(indszDist(sz,:)','Kernel','Kernel','epanechnikov');

       corrXvals = -1:0.01:1;
       yDist{sz,d} = pdf(probDen{sz,d} ,corrXvals);
       yCumdist = cdf(probDen{sz,d},corrXvals);
       idxNull = find(yCumdist >= 0.95);
       if(isempty(idxNull))
            corrNull{sz,d} = 0;   %all weights were most likely set to zero
       else
            corrNull{sz,d} = corrXvals(idxNull(1));
       end
       figNum = 100 + (1-d)*notElim_sz + sz;
       figure(100+sz)
       plot(corrXvals, yDist{sz,d});
%        hist(indszDist(sz,:),20)
       hold on;
%        vline(corrNull{sz,d},'r:',[num2str(lambda(d))])
%only do this stuff for the first lambda value
        if d == numLam
%            hline(0.95,'r:',['Null Corr = ' num2str(corrNull)])
           title(['Fit Null Model Test Seizure No.' num2str(sz) ' For Various Shrinkage Parmater \lambda'])
           vline(testCorr(sz),'b-',['Test Corr = ' num2str(testCorr(sz))])
           legend(lamdaLegend,'location','northwest')
           plotName = ['H:\jaredwil\Lasso Results\permutationTest\' pt{i} '\nullCorrDist_sz' num2str(sz)];
           saveas(gcf,plotName,'jpg')
        end

   end
  
end



%create a struct to save the results
%perm Corr and lambda are index by lambda value
%pdObject,pdf,and nullCorr are all index by sz and lambda values
% permInfo = struct('lambda',lambda,'pdObject',probDen,'pdf',yDist,'nullCorr',corrNull);

permInfo = struct('permCorr',permCorr2,'lambda',lambda,'pdObject',probDen,'pdf',yDist,'nullCorr',corrNull);
infoLabel = ['H:\jaredwil\Lasso Results\permutationTest\' pt{i} '\permInfo.mat'];
save(infoLabel,'permInfo')

%nul distrubution of correlation for ALL SZ 
%this is the correlation distribution for each lambda value
% corrDist = cell(numLam,1);
% for d = 1:numLam
%    xTmp = permCorr2(:,:,d) ;
%    corrDist{d} = xTmp(:);
%    
%    %fit a normal distribution
%    pd = fitdist(xTmp(:),'Normal');
%    corrXvals = -1:0.01:1;
%    yDist = pdf(pd,corrXvals);
%    yCumdist = cdf(pd,corrXvals);
%    
% 
%    
%    figure(d+40)
%    plot(corrXvals,yDist);
%    hold on;
%    plot(corrXvals,yCumdist)
%    hline(0.95,'r:',['Null Corr = ' num2str(corrNull)])
%    title(['Fit Null Model For \lambda = ' num2str(lambda(d))])
% 
% %     [C1, C1x] = hist(xTmp(:),100);
% %     C1 = C1./sum(C1);
% %     plot(C1x,C1)
% %     axis('tight')
% %     title(['Sz Correlation Null \lambda = ' num2str(lambda(d))])
% %     hold on;
% %     [C2, C2x] = hist((f~= 0),10);
% %     C2 = C2./sum(C2);
% %     bar(C2x,C2)
%    
% end
    
end