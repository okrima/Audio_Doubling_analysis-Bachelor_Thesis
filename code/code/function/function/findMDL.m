function [ MDL, errore ] = findMDL( signal, maxOrder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    errore = zeros(size(signal));

    MDL = zeros(maxOrder, size(signal,2));

    for ordine = 1 : maxOrder
        a = lpc(signal, ordine);
        for c = 1 : size(signal,2) - 1
            errore(:,c) = filter(a(c,:), 1, signal(:,c));
        end
        MDL(ordine,:) = MDLtesting(errore,ordine);
    end
    
end
