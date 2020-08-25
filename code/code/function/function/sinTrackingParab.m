function [ Out ] = sinTrackingParab( sr, nfft, locP, offset )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    Out = zeros(1, size(locP, 2));
    tick = (sr/nfft);
    countRiga = 0;
    slope = 0.0;
    for c = 1 : size(locP, 2)
        r = 1;
        q = c;
        listInd = find(locP(:,c));   %Creo lista Picchi per frame
        if isempty(listInd)
            continue;
        end
        while r <= length(listInd)
            n_frame = max(round((offset/tick)) + round(listInd(r) * slope), 1);
            c = q;
            % If is Garbage Peak
            if all(locP(listInd(r) : min(listInd(r) + n_frame, size(locP, 1)), min(c + 1, size(locP, 2))) == 0) && ...
                    all(locP(max(listInd(r) - n_frame, 1) : listInd(r), min(c + 1, size(locP, 2))) == 0)
                locP(listInd(r), c) = 0;
                r = r + 1;
                continue;                
            else
                countRiga = countRiga + 1;
                while true && c <= size(locP, 2)
                    n_frame = max(round((offset/tick)) + round(listInd(r) * slope), 1);
                    arrayP = max(find(locP(max(listInd(r) - n_frame , 1) : listInd(r), min(c + 1, size(locP, 2))))); %Find Previous bin Peak
                    arrayP = n_frame - arrayP + 1;
                    if isempty(arrayP)
                        arrayP = +inf;
                    end
                    arrayN = find((locP(listInd(r) : min(listInd(r) + n_frame, size(locP, 1)), min(c + 1, size(locP, 2)))), 1);
                    arrayN = arrayN - 1;
                    if isempty(arrayN)
                        arrayN = +inf;
                    end
                    passo = min(arrayN, arrayP);
                    %If find Fork (bivio equidistante - possibile biforcazione)
                    if arrayN == arrayP && passo ~= 0 && passo ~= inf
                        countP = 0;
                        countN = 0;
                        tmpPasso = passo;
                        serpe = max(listInd(r) - tmpPasso, 1);
                        countP = countP + 1;
                        tmpC = c+1;
                        %%Itero biforcazione precedente
                        while true && tmpC <= size(locP, 2)
                            n_frame = max(round((offset/tick)) + round(serpe * slope), 1); 
                            %Find Previous bin Peak
                            tmpPrev = max(find(locP(max(serpe - n_frame , 1) : serpe, min(tmpC + 1, size(locP, 2))))); 
                            tmpPrev = n_frame - tmpPrev + 1;
                            if isempty(tmpPrev)
                                tmpPrev = +inf;
                            end
                            %Find Successive bin Peak
                            tmpNext  = find((locP(serpe : min(serpe + n_frame, size(locP, 1)), min( tmpC+ 1, size(locP, 2)))), 1);
                            tmpNext  = tmpNext  - 1;
                            if isempty(tmpNext )
                                tmpNext  = +inf;
                            end
                            tmpPasso = min(tmpNext , tmpPrev);
                            if tmpPasso == inf
                                break;
                            end
                            if tmpPasso == tmpPrev == tmpNext
                                break;                               
                            elseif tmpPasso == tmpPrev
                                serpe = max(serpe - tmpPasso, 1);
                            elseif tmpPasso == tmpNext
                                serpe = max(serpe + tmpPasso, 1);
                            end
                            countP = countP + 1;
                            tmpC = tmpC+1;
                            if tmpC == size(locP, 2) +1
                                break;
                            end
                        end
                        tmpPasso = passo;
                        serpe = max(listInd(r) + tmpPasso, 1);
                        countN = countN + 1;
                        tmpC = c+1;
                        while true && tmpC <= size(locP, 2)
                            n_frame = max(round((offset/tick)) + round(serpe * slope), 1); 
                            %Find Previous bin Peak
                            tmpPrev = max(find(locP(max(serpe - n_frame , 1) : serpe, min(tmpC + 1, size(locP, 2))))); 
                            tmpPrev = n_frame - tmpPrev + 1;
                            if isempty(tmpPrev)
                                tmpPrev = +inf;
                            end
                            %Find Successive bin Peak
                            tmpNext = find((locP(serpe : min(serpe + n_frame, size(locP, 1)), min(tmpC + 1, size(locP, 2)))), 1);
                            tmpNext = tmpNext - 1;
                            if isempty(tmpNext)
                                tmpNext = +inf;
                            end
                            tmpPasso = min(tmpNext, tmpPrev);
                            if tmpPasso == inf
                                break;
                            end
                            if tmpPasso == tmpPrev == tmpNext
                                break;
                            elseif tmpPasso == tmpPrev
                                serpe = max(serpe - tmpPasso, 1);
                            elseif tmpPasso == tmpNext
                                serpe = max(serpe + tmpPasso, 1);
                            end
                            countN = countN + 1;
                            tmpC = tmpC+1;
                            if tmpC == size(locP, 2)
                                break;
                            end
                        end
                        if countP > countN
                            Out(countRiga, c) = (locP(listInd(r), c) - 1) * sr/nfft;
                            locP(listInd(r), c) = 0;
                             listInd(r) = max(listInd(r) - passo, 1);                              
                        elseif countP < countN
                            Out(countRiga, c) = (locP(listInd(r), c) - 1) * sr/nfft;
                            locP(listInd(r), c) = 0;
                            listInd(r) = min(listInd(r) + passo, size(locP, 1)); 
                        else
                            Out(countRiga, c) = (locP(listInd(r), c) - 1) * sr/nfft;
                            locP(listInd(r), c) = 0;
                            listInd(r) = max(listInd(r) - passo, 1); 
                        end
                    else
                        if passo == 0
                            Out(countRiga, c) = (locP(listInd(r), c) - 1)*sr/nfft;
                            locP(listInd(r), c) = 0;
                        elseif passo == arrayP && passo ~= inf %&& locP(max(listInd(r) - passo, 1), max(c - 1, 1)) == 0
                            Out(countRiga, c) = (locP(listInd(r), c) - 1)*sr/nfft;
                            locP(listInd(r), c) = 0;
                            listInd(r) = max(listInd(r) - passo, 1);
                        elseif passo == arrayN && passo ~= inf %&& locP(min(listInd(r) + passo, size(locP, 1)), max(c - 1, 1)) == 0
                            Out(countRiga, c) = (locP(listInd(r), c) - 1)*sr/nfft;
                            locP(listInd(r), c) = 0;
                            listInd(r) = min(listInd(r) + passo, size(locP, 1));
                        elseif passo == inf
                            Out(countRiga, c) = (locP(listInd(r), c) - 1)*sr/nfft;
                            locP(listInd(r), c) = 0;
                            break;
%                         else
%                             Out(countRiga, c) = (locP(listInd(r), c) - 1)*sr/nfft;
%                             locP(listInd(r), c) = 0;
                        end 
                    end
                    c = c + 1;
                end 
            end
            r = r + 1;
        end
    end
end

