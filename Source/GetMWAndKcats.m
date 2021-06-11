%much of this code is copied from convertToEnzymeModel
function MWDivKcatsRes = GetMWAndKcatForReaction(uniprots, kcats, swissprot, standardMW)

swissProtHashMap = containers.Map(swissprot(:,1),1:length(swissprot(:,1)));

[nRxns,~] = size(kcats);
MWDivKcatsRes = NaN(nRxns,1);
[~,n]   = size(uniprots);
y       = 0;
N1 = 0;

for i = 1:length(MWDivKcatsRes)

    %loop through all uniprots - not sure this is needed here, can probably
    %remove this later
    for j = 1:n
        if (isempty(uniprots{i,j}))
            break;
        end
        %Update vector enzymes and count the number of isozymes (x):
        if kcats(i,j) > 0 % if kcat=0 no mets will be added
            uniprots{i,j} = strsplit(uniprots{i,j},' ');
        end
    end

    if ~isempty(uniprots{i,1})
        x = 0;
        MWDivKcats = NaN(n, 1);
        for j = 1:n
            if isempty(uniprots{i,j}) %this means we have reached the end of the list
                break;
            end
            if ~isempty(uniprots{i,j}) && kcats(i,j) > 0 % if kcat=0 no mets will be added
                x         = x+1;
                kcat = kcats(i,j);
                protNames = uniprots{i,j};
                %so, we reason that if the MW is not found, use let's a standard value
                %if a reaction requires 3 proteins, and we only find MWs for 2
                %of them, a standard value is a better guess than zero
                MWs = repmat(standardMW, length(protNames),1);
                %MWs = NaN(length(protNames),1);%this is temporary to compare with the other EC model
                %find MW for the proteins
                %this is likely super slow...

                indices = values(swissProtHashMap,protNames); %so, all must be found; I think they should be found? If not, replace with the loop below
                for k = 1:length(protNames)
                    if swissprot{indices{k},5} ~= 0
                        MWs(k) = swissprot{indices{k},5}/1000;	%g/mmol
                    end
    %                for m = 1:length(swissprot)
    %                    %Molecular Weight:
    %                    if strcmp(protNames{k},swissprot{m,1}) && swissprot{m,5} ~= 0
    %                        MWs(k) = swissprot{m,5}/1000;	%g/mmol
    %                        break;
    %                    end
    %                end
                end

                %these checks are temporary, just remove later
                MWs(isnan(MWs)) = [];
                if (~isempty(MWs))
                    MWDivKcats(j) = sum(MWs/kcat);
                end
            end
        end
        if (x > 0)
            if (any(~isnan(MWDivKcats(1:x))))
                MWDivKcatsRes(i) = min(MWDivKcats(~isnan(MWDivKcats)));
            end
        end
    end
    if rem(i,100) == 0 || i == nRxns 
        fprintf('.')
    end

end
end