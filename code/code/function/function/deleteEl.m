function [ Y ] = deleteEl( X, F )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    Y = zeros(size(X, 1), size(X, 2));
    for track = 1 : size(X, 1)
        for frame = 1 : size(X, 2)
            if X(track, frame) < F
                continue;
            else
                Y(track, frame) = X(track, frame);
            end
        end
    end

end

