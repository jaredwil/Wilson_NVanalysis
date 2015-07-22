function [ feats ] = szPred_winFeatExt( y , fs)
% [ feats ] = szPred_winFeatExt( y )
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
deltaB = [0 4];


%%Begin Function
    [PSD,F]  = pwelch(y,ones(length(y),1),0,length(y),fs,'psd');

    Amp = AmpFn(y);
    LL = LLFn(y);
    nlEng = EnergyNl(y);
    alpha = bandpower(PSD,F,alphaB,'psd');
    beta  = bandpower(PSD,F,betaB,'psd');
    delta = bandpower(PSD,F,deltaB,'psd');
    gamma = bandpower(PSD,F,gammaB,'psd');
    theta = bandpower(PSD,F,thetaB,'psd');

    feats = [Amp LL nlEng alpha beta delta gamma theta];

end

