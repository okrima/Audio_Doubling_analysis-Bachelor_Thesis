function [ Freq ] = diffPhase( Xph, nfft, noverlap, sr )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    freqPerBin = sr / nfft;
    
    Freq = zeros(size(Xph, 1), size(Xph, 2));     

    for bin = 1 : size(Xph, 1)       
%         lastPhase = Xph(bin, 1);
        lastPhase = Xph(bin, 1);
        expect = 2 * pi* (bin - 1) / noverlap; 
        expect = ModSignedPi(expect);
        for frame = 1 : size(Xph, 2) - 1              
            % Differenza di Fase
            dPhase = Xph(bin, frame + 1) - lastPhase;
            dPhase = ModSignedPi(dPhase);
            lastPhase = Xph(bin, frame + 1);
            
            % Tolgo quella che ci si aspetta
            DdPhase = dPhase - (expect);
            
            DPhasePi = ModSignedPi(DdPhase); 

            dBin =  DPhasePi * noverlap / (2*pi);   
            
            Freq(bin, frame+1) = (bin - 1) * freqPerBin + (dBin * freqPerBin);                      
                                
        end
    end
end

