function [ Mid, Side ] = ms( L, R )
% ms : Convert LR in MS 
%   Detailed explanation goes here

    q = sqrt(2);
    Mid(:,1) = (L + R) / q;
    Side(:,1) = (L - R) / q;

end

