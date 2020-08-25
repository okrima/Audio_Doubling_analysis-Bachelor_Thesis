function [ loc, val ] = interParab( Xinput, Rif )
%interParab : Interpolazione Parabolica
%   Detailed explanation goes here
%   Interpolazione parabolica sul bin di valore massimo (peak interolation)
%   Xinput, Rif, Riffase -> Elementi dei picchi che superano la MAF(FFT),
%   Riferimento alla magnitudine, Riferimento alla fase
%   loc, val, valph -> Loczione del vertice, Valore del vertice,
%   Valore della relativa fase
        loc = NaN(size(Xinput,1), size(Xinput,2));
        val = NaN(size(Xinput,1), size(Xinput,2));
        valph = NaN(size(Xinput,1), size(Xinput,2));
        for c = 1:size(Xinput,2)
            for r = 1:size(Xinput,1)
                if isnan(Xinput(r,c))
                   continue; 
                end
                if (r == size(Xinput,1) || r == 1)
                    loc(r,c) = Xinput(r,c);
                    val(r,c) = Rif(Xinput(r,c),c);
                    continue;  
                elseif Xinput(r,c) > 1
                    centroPos = Xinput(r,c);
                    precPos = centroPos - 1;
                    succPos = centroPos + 1;
                    centroVal = Rif(centroPos,c);
                    precVal = Rif(precPos,c);
                    succVal = Rif(succPos,c);
                    vertice = centroPos + (0.5*((precVal - succVal)/(precVal - 2*centroVal + succVal)));
                    loc(r,c) = vertice;
                    ampiezza = centroVal - 0.25*((vertice - centroPos)*(precVal - succVal));
                    val(r,c) = ampiezza; 
%                     valph(r,c) = interp1([1 3],[Riffase(precPos,c) Riffase(succPos,c)],1.5);               
                end
            end
        end
end

