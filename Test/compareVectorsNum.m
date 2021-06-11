function compareVectorsNum(v1,v2, str)
if size(v1) == size(v2)
    %handle NaNs
    sel = isnan(v1);
    sel2 = isnan(v2);
    nanOk = 1;
    if sum(sum(sel)) > 0
        nanOk = all(all(sel == sel2));
    end
    if nanOk 
        v1(sel) = 0; %just get rid of NaNs, they cause problems
        v2(sel2) = 0;
        if all(all(v1 == v2))
            disp([str ': ok'])   
        else
            disp([str ': values differ'])
        end
    else
        disp([str ': NaNs differ'])
    end
else
    disp([str ': size differ'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%