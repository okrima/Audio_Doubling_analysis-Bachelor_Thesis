function [ Y ] = ModSignedPi( DdPhase )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Q = fix(DdPhase / pi);     
            
    if mod(Q,2) ~= 0
        if Q >= 0
            Q = Q + 1;
        else
            Q = Q - 1;
        end
    end

Y = DdPhase - Q * pi; 

end

