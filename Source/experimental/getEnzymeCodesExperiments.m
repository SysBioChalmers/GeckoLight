%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_data = getEnzymeCodes(model,action)
% Retrieves the enzyme codes for each of the reactions for a given genome
% scale model (GEM).
%
% INPUT:    *model:      a GEM file (.mat format)
%           *action:     Response action if multiple proteins with
%                        different EC numbers are found for a given gene in
%                        a metabolic reaction.
%                           - 'display' Displays all found multiplicities.
%                           - 'ignore'  Ignore multiplicities and use the
%                              protein with the lowest index in the database.
%                           - 'add'     Adds all the multiple proteins as
%                                       isoenzymes for the given reaction.
%           
% OUTPUTS:  model_data, which contains:
%           *model:      The standardized GEM
%           *substrates: Substrates associated for each rxn
%           *products:   Products associated, when rxn is reversible
%           *uniprots:   All possible uniprot codes, for each rxn
%           *EC_numbers: All possible EC numbers, for each uniprot
%           *count(1):   #rxns with data from Swissprot
%           *count(2):   #rxns with data from KEGG
%           *count(3):   #exchange/transport rxns with no GPRs
%           *count(4):   #other rxns
% 
% Benjamin Sanchez. Last edited: 2017-03-05
% Ivan Domenzain.   Last edited: 2018-09-07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model_data = getEnzymeCodes(model,action)

if nargin<2
    action = 'display';
end

fprintf('Retrieving EC numbers...')

%Standardize grRules to avoid wrong enzyme codes assignments to reactions
[grRules,~]   = standardizeGrRules(model,true);
model.grRules = grRules;

data      = load('../../databases/ProtDatabase.mat');
swissprot = data.swissprot;
kegg      = data.kegg;

swissprot = standardizeDatabase(swissprot);
kegg      = standardizeDatabase(kegg);

DBprotSwissprot     = swissprot(:,1);
DBgenesSwissprot    = flattenCell(swissprot(:,3));
DBecNumSwissprot    = swissprot(:,4);
DBMWSwissprot       = swissprot(:,5);

DBprotKEGG          = kegg(:,1);
DBgenesKEGG         = flattenCell(kegg(:,3));
DBecNumKEGG         = kegg(:,4);
DBMWKEGG            = kegg(:,5);

[m,n]      = size(model.S);
substrates = cell(n,20);
products   = cell(n,20);
uniprots   = cell(n,20);
Genes      = cell(n,20);
EC_numbers = cell(n,20);
MWs        = zeros(n,20);
isrev      = zeros(n,1);
count      = zeros(4,1);
rgmat      = full(model.rxnGeneMat);
conflicts  = cell(1,4);
tic
for i = 1:n
    ks  = 1;
    kp  = 1;
    dir = 0;
    inv = 0;
    DB = '';
    %Save the substrates and products (if rxn is reversible):
    for j = 1:m
        if model.S(j,i) < 0 && model.ub(i) > 0
            substrates{i,ks} = model.metNames{j};
            ks  = ks+1;
            dir = 1;
        elseif model.S(j,i) > 0 && model.lb(i) < 0
            products{i,kp} = model.metNames{j};
            kp  = kp+1;
            inv = 1;
        end
    end
    if rem(i,100) == 0 || i == n
        fprintf('.')
    end
end
toc %Elapsed time is 639.305883 seconds.
tic
for i = 1:n
    ks  = 1;
    kp  = 1;
    dir = 0;
    inv = 0;
    DB = '';
    %Save the substrates and products (if rxn is reversible):
    %This code was opimized (a loop through every element was replaced)
    %It was tested and produces the same results, but roughly 2000 times faster
    if model.ub(i) > 0
       sbstrs = model.metNames(model.S(:,i) < 0); 
       substrates(i,1:length(sbstrs)) = sbstrs;
    elseif model.lb(i) < 0 
       prds = model.metNames(model.S(:,i) > 0); 
       products(i,1:length(prds)) = prds;
    end
    if rem(i,100) == 0 || i == n
        fprintf('.')
    end
end
toc %Elapsed time is 0.293855 seconds, so roughly 2000 times faster.

%Test to build an index from gene to prot
tic
[x,y] = size(DBgenesSwissprot);
genesForIndex = reshape(DBgenesSwissprot, x*y, 1);
genesForIndex = genesForIndex(~cellfun(@isempty, genesForIndex));
genesForIndex = unique(genesForIndex); %18360
%genesForIndexCat = categorical(genesForIndex,genesForIndex);
geneIndex = cell(length(genesForIndex),1);
geneHashMap = containers.Map(genesForIndex,1:length(genesForIndex));
protIndices = 1:length(DBgenesSwissprot(:,1));
for i = 1:10
    tmp1 = DBgenesSwissprot(:,i);
    sel = ~cellfun(@isempty, tmp1);
    indices = cell2mat(values(geneHashMap,tmp1(sel)));
    protIndicesSel = protIndices(sel);
    for j = 1:length(indices)
        geneIndex{indices(j)} = [geneIndex{indices(j)};protIndicesSel(j)]; 
    end
end
toc %Elapsed time is 0.391342 seconds, ok.


tic
for i = 1:n
    if ~isempty(model.grRules{i})
        %Find match in Swissprot:
        [new_uni,new_EC,new_MW,newGene,multGenes] = findInDBOpt(model.grRules{i},DBprotSwissprot,DBgenesSwissprot,DBecNumSwissprot,DBMWSwissprot,geneIndex,geneHashMap);
    end
    if rem(i,100) == 0 || i == n
        fprintf('.')
    end
end
toc

