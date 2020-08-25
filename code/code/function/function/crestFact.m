function [ Xcrest ] = crestFact( Xpicco, Xrms )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%   Calcolo fattore di cresta
%   Xpicco, Xrms -> Valore picco e rms
%   Xcrest -> Valore di crest factor

    Xcrest = Xpicco-Xrms;
    
end

