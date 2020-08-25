function [ F ] = Cepstrum( X )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    F = zeros(size(X, 1), size(X, 2));

    for frame = 1 : size(X, 2)
        C = X(:,frame);
        Pow = log10(abs(C).^2);
        F(:,frame) = ifft((Pow)); 
    end

end
   
