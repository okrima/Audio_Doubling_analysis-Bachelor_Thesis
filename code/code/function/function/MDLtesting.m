function [ x ] = MDLtesting( errore, ord )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    N = size(errore,1);     % number of data points fitted
    V = mean(errore.^2);     % loss function;
    x = V * (1 + ord * log(N) / N);     % Rissanen's Minimum Description Length


end

