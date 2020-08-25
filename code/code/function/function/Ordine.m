function [ S ] = Ordine(Sig, maxOrder)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    [ MDL, ~ ] = findMDL( Sig, maxOrder);

    [ ~, S ] = min( MDL );
    
    S = round(mean(S));
    
end

