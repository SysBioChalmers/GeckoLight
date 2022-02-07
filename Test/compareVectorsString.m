function compareVectorsString(v1,v2, str)
if size(v1) == size(v2)
    if all(all(strcmp(v1,v2)))
        disp([str ': ok'])   
    else
        disp([str ': string values differ'])
    end
else
    disp([str ': size differ'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%