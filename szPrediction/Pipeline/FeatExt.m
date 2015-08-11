function [ feats ] = FeatExt( y , fs)
% [ feats ] = FeatExt( y )
% Input to this function is a time window y and sampling freq (fs). Output
% is array of features detailed below
% FEATURES TO BE EXTRACTED - This will be updated as algorithm progresses
%           - Area
%           - Line Length
%           - Non-Linear (Teager) Energy
%           - Spectral Power:
%                       -- Alpha
%                       -- Beta
%                       -- Delta
%                       -- Gamma
%                       -- Theta

%%
% Anonymous Feature Functions
AmpFn = @(x) mean(abs(x),1);
AreaFn = @(x) sum(abs(x),1);
EnergyFn = @(x) mean(x.^2,1);
ZCFn = @(x) sum((x(1:end-1,:)>repmat(mean(x),size(x,1)-1,1)) & x(2:end,:)<repmat(mean(x),size(x,1)-1,1) | (x(1:end-1,:)<repmat(mean(x),size(x,1)-1,1) & x(2:end,:)>repmat(mean(x),size(x,1)-1,1)));
LLFn = @(x) mean(abs(diff(x)));
LLFn2 = @(X, winLen) conv2(abs(diff(X,1)),  repmat(1/winLen,winLen,1),'same');

%mean windowed nonlinear energy -- teager energy. 
EnergyNl = @(x) mean(bsxfun(@minus,x(2:end-1,:).^2, x(1:end-2,:).*x(3:end,:)),1);

%%
%Define Constants 
%spectral bands in Hz
alphaB = [8 12];
betaB = [12 27];
gammaB = [27 45];
thetaB = [3 8];
deltaB = [0.1 4];

% for upper triangle indexing
ix = triu( ones(size(y,2)), 1 );

%define options for eigs function 
opts.tol = 1e-3;  %change tolerance

%%Begin Function
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Spectral Featrues
%  tic;
    [PSD,F]  = pwelch(y,ones(length(y),1),0,length(y),fs,'psd'); 
    PSD(F > 47,:) = [];
    F(F > 47) = [];
    alpha = bandpower(PSD,F,alphaB,'psd');
    beta  = bandpower(PSD,F,betaB,'psd');
    delta = bandpower(PSD,F,deltaB,'psd');
    gamma = bandpower(PSD,F,gammaB,'psd');
    theta = bandpower(PSD,F,thetaB,'psd');
% toc;

%average across freq bins
    PSDnorm = bsxfun(@rdivide,bsxfun(@minus,PSD,mean(PSD,2)), std(PSD,[],2) + 2e-13);
    %nonlinear
    %cross correlation between channels
    psdCorr = corrcoef(PSDnorm(2:end,:));  %omit zero hurtz
%     psdEigs = eigs(psdCorr,opts)';              %only return the six highest
    
    % flatten upper triangle of correlation matrix
    psdCorr = triu( psdCorr, 1 );
    psdCorr = psdCorr( find(ix) )';
%     save('fuckedupWindow.mat','y');  %used to find the window that caused
%                                      %the error
    assert( sum(isnan(psdCorr)) == 0 , 'Nans');

%is computing fft PSD faster??? --NO     
%     tic;
%     Pxx = abs(fft(y)).^2;
%     pFreq = linspace(0,fs,length(y)/2 + 1);
%     
%     alpha = bandpower(Pxx(1:length(y)/2 + 1,:),pFreq,alphaB,'psd');
%     beta  = bandpower(Pxx(1:length(y)/2 + 1,:),pFreq,betaB,'psd');
%     delta = bandpower(Pxx(1:length(y)/2 + 1,:),pFreq,deltaB,'psd');
%     gamma = bandpower(Pxx(1:length(y)/2 + 1,:),pFreq,gammaB,'psd');
%     theta = bandpower(Pxx(1:length(y)/2 + 1,:),pFreq,thetaB,'psd');
%     toc;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %Temporal Features
    %amplitude
    Amp = AmpFn(y);
    %linelength
    LL = LLFn(y);
    %teager energy
    nlEng = EnergyNl(y);
    

    %skewness and kurtosis
    fxskew = skewness(y);
    fxkurt = kurtosis(y);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%nonlinear
    %%cross correlation between channels
    ycorr = corrcoef(y);
%     yEigs = eigs(ycorr,opts)';  %only return the six highest
    %flatten correlation matrix
    ycorr = triu( ycorr, 1 );
    fycorr = ycorr( find(ix) )';
    
    
    %%channel covariance 
    ycov = cov(y);
    % flatten upper triangle of covariance matrix
    ycov = triu( ycov, 1 );
    fycov = ycov( find(ix) )';
    
    %KFD fractal Dim
    
    %Hurst exponenst
%     KFD = zeros(1,size(y,2));
    HFD = zeros(1,size(y,2));  %higuchi fractal dimension
    kmax = 2; %kmax = 2
    for ch = 1:size(y,2)
%         KFD(ch) = Katz_FD(y(:,ch));
        HFD(ch) = Higuchi_FD(y(:,ch),kmax);   
    end
    
    feats = [Amp LL nlEng alpha beta delta gamma theta fxskew fxkurt psdCorr fycorr fycov HFD];

end

