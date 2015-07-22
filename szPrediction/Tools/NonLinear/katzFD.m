function [ D ] = katzFD( x )
%UNTITLED Summary of this function goes here
% By: Jared Wilson
% 7/14/2015 
%   This function takes in a signal and calculates the fractal dimension
%   using the Katz algorithm. 
%Inputs:
% x - signal to be analyzed
%outputs:
% D - Fractal Dimension calculated using Katz FD algorithm
%%
%annonomous funcitons needed
avgDif = @(x) mean(sqrt(1+diff(x).^2));
curveLen = @(x) sum(sqrt(1+diff(x).^2));  %use 1 becuase the difference will 
                                        %always be one between consecutive numbers
maxDist = @(x) max(sqrt( bsxfun(@plus, reshape(((2:length(x)) - 1).^2,size(x(2:end))), ...
    bsxfun(@minus,repmat(x(1),size(x(2:end))),x(2:end)).^2 )));

L = curveLen(x);
a = avgDif(x);
d = maxDist(x);


n = L/a;

D = log10(n)/(log10(d/L) + log10(n));
%%start
end

