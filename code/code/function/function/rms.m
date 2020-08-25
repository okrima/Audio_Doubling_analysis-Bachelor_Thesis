function [ Xrms ] = rms( X )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%   X -> Input signal
%   Xrms - Output, Valore rms
    Xrms = sqrt(mean(X.^2,1));

end

