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

%%
%begin function
%loop through all pts
for i = 1:numel(pt);  

    
    featLab = ['C:\Users\Jared\Dropbox\NVanalysis_data\SzPred_data\Res LassoLib\' pt{i} '_bestLasso5.mat'];

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
    
    featInfo_all(i,:) = numNotZero;

    %is there a way to find out which channels contain important information
    %for regression????
    %For ch-by-ch features this is very easy... for the cov and corr matrices
    %this is a bit more challenging.
    rowCh = [];
    colCh = [];
    for i = 1:15
        rowCh = [rowCh 1:i];
        colCh = [colCh ones(1,i)*(i+1)];
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


end


%plot some stuff
plotnumeric(featInfo_all);