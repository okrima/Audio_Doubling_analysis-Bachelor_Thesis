function [ Y ] = zero2NaN( X )
%zero2NaN Convert element 0 to NaN
%   Detailed explanation goes here

Y = NaN(size(X, 1), size(X, 2));

    for c = 1 : size(X, 2)
        for r = 1 : size(X, 1)
            if X(r,c) == 0
                continue;
            else
                Y(r,c) = X(r,c);
            end
        end
    end

end

