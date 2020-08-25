function [ Y ] = PFD( Mag, Thres )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Normalizzazione magnitudine tra 0 e 1
    Mag = amp2db(Mag, -96) - 96;
    Mag = Mag./max(Mag(:));

    % Detection
    Det = 20*log10(Mag(:,1:end - 1)./(Mag(:,2:end))) > Thres;

    % Percussive feature
    Y = sum(Det, 1);
    
end
