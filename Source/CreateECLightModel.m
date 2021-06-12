function ecLightModel = CreateECLightModel(model, fillInMissingGPRs, speciesAdapter)
 
%There are four major points with this function:
%1. reduce the size of the ec model, since that is a huge problem when you
%try to simulate several cell types together - it is just too slow.
%2. There are a lot of missing GPRs for reactions - here we fill in a
%standard cost for such reactions.
%3. It is much faster - 1-2 minutes instead of hours.
%4. The standard ec model gives the solver some numerical problems - it 
%sometimes fails to solve the FBA. The light model seems to work fine.

%much of this code is copied from various parts of Gecko and ecModels, 
%such as the function enhanceGEM

if nargin < 3
	speciesAdapter = HumanGEMAdapter();
end
if nargin < 2
	fillInMissingGPRs = true;
end


%Get model-specific parameters
parameters = speciesAdapter.getParameters();

%Remove blocked rxns + correct model.rev:
[model,name,modelVer] = preprocessModelLt(model,'ecLight','1');

fprintf('\n==================')
fprintf('\nGenerating light ecModel:')
fprintf('\n==================\n')

model_data = getEnzymeCodesOpt(model, speciesAdapter); %Elapsed time is 38.600622 seconds, I think this could be acceptable.
fprintf('\n')


%Load BRENDA data:
[KCATcell, SAcell] = loadBRENDAdataLt(speciesAdapter);
kcats      = matchKcatsOpt(model_data,parameters.org_name, KCATcell, SAcell, speciesAdapter);
fprintf('\n')
%code here copied from readKcatData

%Get kcat value for both directions:
Fkcat = kcats.forw.kcats;
Bkcat = kcats.back.kcats;
rev   = logical(model_data.model.rev);
kcats = [Fkcat;Bkcat(rev,:)];

%Update uniprots with both directions:
uniprots = [model_data.uniprots; model_data.uniprots(rev,:)];

%Update matched genes with both directions:
matchedGenes = [model_data.matchedGenes; model_data.matchedGenes(rev,:)];

%Convert to irreversible model with RAVEN function (will split in 2 any reversible rxn):
model = convertToIrrev(model_data.model);


%What we want to find for each reaction
%is the min value of MW/kcat when comparing all parallel reactions (i.e. OR in the GR rules)
%The plan is then to add one reaction with prot_pool_exchange, one
%metabolite called prot_pool, and a stochiometric coefficient of MW/kcat of
%that for each reaction. For complexes, the MW are summed up, since the kcat is the same.

%Load databases:
data      = load(speciesAdapter.getFilePath('ProtDatabase.mat')); %this is loaded twice, (also in a function above), a bit unnecessary
swissprot = data.swissprot;

standardMW = median([swissprot{:,5}])/1000;%the median of all the proteins

%profile on -history
%profile off
%p = profile('info')
%[srt,I] = sort([p.FunctionTable.TotalTime],2, 'ascend');
%fnksTmp = {p.FunctionTable.FunctionName};
%fnks = fnksTmp(I);
%fnks{54}
%srt(54)

fprintf('Creating Gecko light model')

MWDivKcats = GetMWAndKcats(uniprots,kcats, swissprot, standardMW);
fprintf('\n');

%so, what we do now is to add the reaction prot_pool_exchange and the metabolite prot_pool
rxnsToAdd = struct();
rxnsToAdd.rxns = {'prot_pool_exchange'};
rxnsToAdd.equations = {'=> prot_pool[c]'};
rxnsToAdd.rxnNames = {'prot_pool_exchange'};
metsToAdd = struct();
metsToAdd.mets = {'prot_pool'};
metsToAdd.metNames = {'prot_pool'};
metsToAdd.compartments = {'c'};

model = addMets(model, metsToAdd);

%Now add the protein metabolite to the reactions
metRow = MWDivKcats.';
metRow(isnan(metRow)) = 0;
model.S(length(model.mets),:) = -metRow;

model = addRxns(model, rxnsToAdd, 3);

%set the protein pool constraint
model.ub(strcmp(model.rxns, 'prot_pool_exchange')) = 0.05900217; %This was estimated in the TME modeling project.



if fillInMissingGPRs
    %Now comes the task of filling in a standard protein cost for reactions
    %with missing GPRs
    standardRxnProtCost = median(MWDivKcats(~isnan(MWDivKcats)));%1.2894e-04
    MWDivKcatsFilledIn = [MWDivKcats;1]; %add prot_pool_exchange

    %get info about spontanoeus rxns
    [spont,spontRxnNames] = speciesAdapter.getSpontaneousReactions(model);


    [~,exchangeRxnsIndexes] = getExchangeRxns(model);

    numToFix = 0;
    protPoolIndex = find(strcmp(model.mets, 'prot_pool'));

    for i = 1:(length(MWDivKcats)-1)%avoid the protein_pool_exchange
       if isnan(MWDivKcats(i))%only look at reactions without protein cost
           %Step 1: Skip exchange reactions
           if ismember (i,exchangeRxnsIndexes)
               %disp(strcat(num2str(i),': skipping exchange: ', model.rxns(i), ': ', constructEquations(model, model.rxns(i))));
               continue;
           end
           %Step 2: Skip transport reactions (with missing GPRs)
           %We really only skip the reactions with two metabolites
           %where they are the same but different compartments
           numMets = sum(model.S(:,i) ~= 0);
           if numMets == 2
              mets = model.metComps(model.S(:,i) ~= 0);
              metNames = model.metNames(model.S(:,i) ~= 0);
              if (~strcmp(mets(1), mets(2))) %different compartments
                  if strcmp(metNames{1}, metNames{2}) %the same metabolite, it is a transport reaction
                      %disp(strcat(num2str(i),': skipping transport: ', model.rxns(i), ': ', constructEquations(model, model.rxns(i))));
                      continue;
                  end
              end
           end
           %Step 3: check if this is a spontaneous reaction (last because it is a bit heavier than the other checks)
           rxnName = model.rxns{i};
           if endsWith(rxnName, '_REV')
              rxnName = extractBefore(rxnName, strlength(rxnName)-3);
           end
           spontaneous = spont(strcmp(spontRxnNames, rxnName)); %could potentially be optimized
           if ~isempty(spontaneous)%it shouldn't be
               if ~isnan(spontaneous)%we treat NaN as that it is unknown if the reaction is spontaneous, which is seen as a no
                   if spontaneous
                      %disp(strcat(num2str(i),': skipping spontaneous: ', model.rxns(i), ': ', constructEquations(model, model.rxns(i))));
                      continue; %leave spontaneous reactions without protein cost
                   end
               end
           end
           numToFix = numToFix + 1;
           %disp(strcat(num2str(i),': Assigning std val: ', model.rxns(i), ': ', constructEquations(model, model.rxns(i))));
           %fill in a standard value
           model.S(protPoolIndex, i) = -standardRxnProtCost;
       end
    end
    
    disp(['Filled in ' num2str(numToFix) ' protein usage costs with a standard value due to missing GPRs.'])
end


ecLightModel = model;

end