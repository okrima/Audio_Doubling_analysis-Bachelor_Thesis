function [ c, error, Pgain ] = LinearPredCod( X, ncoef)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    [c, Pgain] = lpc(X, ncoef);
    if isnan(Pgain) 
        c(:) = 0;
        Pgain = 1;
    end
    est_x = filter(c, sqrt(Pgain), X);        % Estimated signal
    error = X - est_x;                        % Residual Signal
  
end

