function [ Y ] = NaN2zero( X )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Y = zeros(size(X, 1), size(X, 2));

    for c = 1 : size(X, 2)
        for r = 1 : size(X, 1)
            if  isnan(X(r,c))
                continue;
            else
                Y(r,c) = X(r,c);
            end
        end
    end

end

