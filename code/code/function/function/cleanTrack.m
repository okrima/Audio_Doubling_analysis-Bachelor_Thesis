function [ locP ] = cleanTrack(locP, t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    r = 1;
    while r < (size(locP, 1) + 1)
        listInd = find(locP(r,:));
        if length(listInd) < t
%             locP(r, listInd(1) : listInd(end)) = 0;
%             Phase(r, listInd(1) : listInd(end)) = 0;
            locP(r,:) = [];
%             Phase(r,:) = [];
            continue;
        end
        r = r + 1;
    end
end

