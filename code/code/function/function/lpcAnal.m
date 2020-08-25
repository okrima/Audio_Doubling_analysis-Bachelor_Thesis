function [ a, g, e, sig] = lpcAnal( s , z, nfft )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    for c = 1 : size(s, 2)
        [ autos , lags ] = xcorr ( s, z) ; 
        autos = autos ( find ( lags == 0) : end ) ; 
        [ a(:,c) , g(:,c) ] = levinson ( autos , nfft) ; 

        e_r(:,c) = filter( a , sqrt ( g ) , z ) ;
        e_l(:,c) = filter( a , sqrt ( g ) , s ) ;
        sig(:,c) = filter( sqrt ( g ) ,a , e ) ; 
    end
 
% 5 autos : Autocorrelation of input signal [ vector ] 
% 6 s : Input signal [ vector ] 
% 7 a : LPC c o e f f i c i e n t s 
% 8 g : Prediction error variance


end

