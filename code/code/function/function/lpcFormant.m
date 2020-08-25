function [ Freq ] = lpcFormant( Coeff, sr )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    r = roots(Coeff);
    r = r(imag(r)>=0);
%   Freq = atan2(imag(r),real(r));
    Freq = angle(r);
    Freq = sort(Freq.*(sr/(2*pi)));   
%     A = (abs(r(I)));
end

