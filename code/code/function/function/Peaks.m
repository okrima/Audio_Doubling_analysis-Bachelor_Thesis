function [ Picchi ] = Peaks( X, soglia, flag )
%UNTITLED9 Summary of this function goes here
%   Peak Detection
%   X -> Input signal
%   soglia -> Soglia di analisi
%   Picchi -> Output
    Picchi = zeros(size(X, 1),size(X, 2));
    if flag == 1
        X = mag2db(X);
        mag = max(X, [], 1);
        mag = max(mag);
        for c = 1 : size(X, 2)
            for r = 1 : size(X, 1)
                if  X(r, c) > mag + soglia
                    Picchi(r, c) = (X(r, c));
                else
                    continue;
                end
            end
        end
    elseif flag == 0
        for c = 1 : size(X, 2)
            for r = 1 : size(X, 1)
                if  (X(r, c)) > soglia
                    Picchi(r, c) = (X(r, c));
                else
                    continue;
                end
            end
        end
    end
end

