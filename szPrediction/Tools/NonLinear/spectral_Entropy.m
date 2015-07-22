function [ entropy ] = spectral_Entropy(y,band,fs)
%[ entropy ] = spectral_Entropy( x )
%   This function calculates the spectral entropy for input signal y for
%   the desired band of interest
% Inputs:
% y  - input time domain signal
% fs - sampling frequency of the signal

% Outputs:
% entropy - calculated entropy
%%%%%%%%%%%%%%  ONE SIDED %%%%%%%%%%%%%

n = length(y);
NFFT = 1024;

yf = fft(y,NFFT,1);
yf = yf(1:NFFT/2,:);
freqBins = fs*linspace(0,0.5,length(yf));

%Remove indicies that are outside the band of interest
bandInds = freqBins >= band(1) & freqBins <= band(2);
yf = yf(bandInds,:);


X = yf.*conj(yf)/(NFFT*n);  %Calculae Power Spectrum Density
Py = bsxfun(@rdivide,X,sum(X + 1e-12,1)); %Create pdf of spectrum

entropy = -sum(Py.*log2(Py+1e-12))/log2(length(Py));

end

