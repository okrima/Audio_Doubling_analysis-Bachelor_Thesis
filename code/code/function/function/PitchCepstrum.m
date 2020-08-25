function [ F0 ] = PitchCepstrum( C, sr )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    F0 = zeros(1, size(C, 2));

    ms1 = floor(sr * 1/250); %  400 Hz
    ms2 = floor(sr * 1/80); %  80 Hz
    
    
    for frame = 1 : size(C, 2)
        Coe = C(:,frame);
        [~ , idx] = max(real(Coe(60:floor(end/2)+1)));
        F0(frame) = sr/(idx + 60);
    end

end

