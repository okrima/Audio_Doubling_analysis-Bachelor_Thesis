function [ Y] = PFD( Mag, Thres )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    % Normalizzazione magnitudine tra 0 e 1
    Mag = amp2db(Mag, -80) -80;
%     Mag = db2mag(Mag);
    Mag = Mag./max(Mag(:));
    % Detection
    Det = mag2db(Mag(:,1:end-1)./Mag(:,2:end)) > Thres;

%     for i = 1 : size(Mag, 1)
%        Det(i,:) = 20*log10(Mag(i, 1:end - 1)./(Mag(i, 2:end))) > Thres;
%     end

    % Percussive feature
    Y = sum(Det, 1);
    
end
