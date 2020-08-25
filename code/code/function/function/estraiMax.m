function [ loc, val ] = estraiMax( Picchi, Z, flag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Estraggo le locazioni del valore massimo tra i bin hanno superato la
%   soglia
%   Xpicchi -> Picchi
%   loc -> Valore della locazione (bin)
    loc = NaN(size(Picchi,1),size(Picchi,2));
    val = NaN(size(Picchi,1),size(Picchi,2));
    if flag == 1
        Z = mag2db(Z);
    end
    for c=1:size(Picchi,2)
        for r=1:size(Picchi,1)
            if Picchi(r,c) ~= 0
                if flag == 1
                    if Z(max(r - 1, 1), c) < Picchi(r,c) && Z(min(r + 1, size(Picchi, 1)), c) < Picchi(r,c)
                            loc(r,c) = r;                    
                            val(r,c) = Picchi(r,c);
                    end
                else
                    if (Z(r, max(c - 1, 1)) < Picchi(r,c) && Z(r, min(c + 1, size(Picchi,2))) <= Picchi(r,c)) || (c==1 && Z(r, min(c + 1, size(Picchi,2))) < Picchi(r,c)) 
                            loc(r,c) = c;                    
                            val(r,c) = Picchi(r,c);
                    end
                end
            end
        end
    end   
end
