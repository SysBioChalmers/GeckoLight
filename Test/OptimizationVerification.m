%%%%%%%%%%%%%%%%%%
%GetEnzymeCodes
%%%%%%%%%%%%%%%%%%

speciesAdapter = HumanGEMAdapter();

tic
model_data = getEnzymeCodesOpt(model, speciesAdapter);
toc %Elapsed time is 38.600622 seconds, I think this could be acceptable.

%check that the result is the same
cd C:\Work\MatlabCode\components\Gecko\GeckoForHuman20210603\GECKO/geckomat/get_enzyme_data
tic
model_data_orig = getEnzymeCodes(model);%Elapsed time is 942.718402 seconds.
toc

length(model_data_orig.substrates) == length(model_data.substrates)%ok
comparisons = false(length(model_data_orig.substrates), 1);
for i = 1:length(model_data_orig.substrates)
    if (isempty(model_data_orig.substrates{i}) & isempty(model_data.substrates{i}))
        comparisons(i) = 1;
    elseif length(model_data_orig.substrates{i}) == length(model_data.substrates{i})
        comparisons(i) = all(strcmp(model_data_orig.substrates{i}, model_data.substrates{i}));
    else
        comparisons(i) = 0;
    end
end
all(comparisons) %ok

length(model_data_orig.products) == length(model_data.products)%ok
comparisons = false(length(model_data_orig.products), 1);
for i = 1:length(model_data_orig.products)
    if (isempty(model_data_orig.products{i}) & isempty(model_data.products{i}))
        comparisons(i) = 1;
    elseif length(model_data_orig.products{i}) == length(model_data.products{i})
        comparisons(i) = all(strcmp(model_data_orig.products{i}, model_data.products{i}));
    else
        comparisons(i) = 0;
    end
end
all(comparisons)  %ok

length(model_data_orig.matchedGenes) == length(model_data.matchedGenes)%ok
comparisons = false(length(model_data_orig.matchedGenes), 1);
for i = 1:length(model_data_orig.matchedGenes)
    if (isempty(model_data_orig.matchedGenes{i}) & isempty(model_data.matchedGenes{i}))
        comparisons(i) = 1;
    elseif length(model_data_orig.matchedGenes{i}) == length(model_data.matchedGenes{i})
        comparisons(i) = all(strcmp(model_data_orig.matchedGenes{i}, model_data.matchedGenes{i}));
    else
        comparisons(i) = 0;
    end
end
all(comparisons)  %ok

length(model_data_orig.uniprots) == length(model_data.uniprots)%ok
comparisons = false(length(model_data_orig.uniprots), 1);
for i = 1:length(model_data_orig.uniprots)
    if (isempty(model_data_orig.uniprots{i}) & isempty(model_data.uniprots{i}))
        comparisons(i) = 1;
    elseif length(model_data_orig.uniprots{i}) == length(model_data.uniprots{i})
        comparisons(i) = all(strcmp(model_data_orig.uniprots{i}, model_data.uniprots{i}));
    else
        comparisons(i) = 0;
    end
end
all(comparisons)  %ok

length(model_data_orig.EC_numbers) == length(model_data.EC_numbers)%ok
comparisons = false(length(model_data_orig.EC_numbers), 1);
for i = 1:length(model_data_orig.EC_numbers)
    if (isempty(model_data_orig.EC_numbers{i}) & isempty(model_data.EC_numbers{i}))
        comparisons(i) = 1;
    elseif length(model_data_orig.EC_numbers{i}) == length(model_data.EC_numbers{i})
        comparisons(i) = all(strcmp(model_data_orig.EC_numbers{i}, model_data.EC_numbers{i}));
    else
        comparisons(i) = 0;
    end
end
all(comparisons)  %ok

length(model_data_orig.MWs) == length(model_data.MWs)%ok
all(all(model_data_orig.MWs == model_data.MWs))%ok

all(model_data_orig.count == model_data.count) %ok

%%%%%%%%%%%%%%%%%%
% matchKcats
%%%%%%%%%%%%%%%%%%


cd C:\Work\MatlabCode\components\Gecko\GeckoForHuman20210603\GECKO/geckomat/get_enzyme_data
tic
kcats_orig      = matchKcats(model_data,parameters.org_name);
toc %Elapsed time is 1521.572380 seconds.

%Load BRENDA data:
[KCATcell, SAcell] = loadBRENDAdata;
cd C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy/GeckoWork
tic
kcats      = matchKcatsOpt(model_data,parameters.org_name, KCATcell, SAcell, 'C:/Work/MatlabCode/components/Gecko/GeckoForHuman20210603/GECKO/Databases/PhylDist.mat');
toc
%compare kcats and kcats_orig:
compareVectorsNum(kcats_orig.forw.kcats, kcats.forw.kcats, 'kcats.forw.kcats');
compareVectorsNum(kcats_orig.forw.org_s, kcats.forw.org_s, 'kcats.forw.org_s');
compareVectorsNum(kcats_orig.forw.rest_s, kcats.forw.rest_s, 'kcats.forw.rest_s');
compareVectorsNum(kcats_orig.forw.org_ns, kcats.forw.org_ns, 'kcats.forw.org_ns');
compareVectorsNum(kcats_orig.forw.rest_ns, kcats.forw.rest_ns, 'kcats.forw.rest_ns');
compareVectorsNum(kcats_orig.forw.org_sa, kcats.forw.org_sa, 'kcats.forw.org_sa');
compareVectorsNum(kcats_orig.forw.rest_sa, kcats.forw.rest_sa, 'kcats.forw.rest_sa');

compareVectorsNum(kcats_orig.back.kcats, kcats.back.kcats, 'kcats.back.kcats');
compareVectorsNum(kcats_orig.back.org_s, kcats.back.org_s, 'kcats.back.org_s');
compareVectorsNum(kcats_orig.back.rest_s, kcats.back.rest_s, 'kcats.back.rest_s');
compareVectorsNum(kcats_orig.back.org_ns, kcats.back.org_ns, 'kcats.back.org_ns');
compareVectorsNum(kcats_orig.back.rest_ns, kcats.back.rest_ns, 'kcats.back.rest_ns');
compareVectorsNum(kcats_orig.back.org_sa, kcats.back.org_sa, 'kcats.back.org_sa');
compareVectorsNum(kcats_orig.back.rest_sa, kcats.back.rest_sa, 'kcats.back.rest_sa');

kcats_orig.tot.queries == kcats.tot.queries
kcats_orig.tot.org_s == kcats.tot.org_s
kcats_orig.tot.rest_s == kcats.tot.rest_s
kcats_orig.tot.org_ns == kcats.tot.org_ns
kcats_orig.tot.rest_ns == kcats.tot.rest_ns
kcats_orig.tot.org_sa == kcats.tot.org_sa
kcats_orig.tot.rest_sa == kcats.tot.rest_sa
kcats_orig.tot.wc0 == kcats.tot.wc0
kcats_orig.tot.wc1 == kcats.tot.wc1
kcats_orig.tot.wc2 == kcats.tot.wc2
kcats_orig.tot.wc3 == kcats.tot.wc3
kcats_orig.tot.wc4 == kcats.tot.wc4
compareVectorsNum(kcats_orig.tot.matrix, kcats.tot.matrix, 'kcats.tot.matrix');

MWDivKcats_orig = NaN(length(model.rxns),1);

tic
for i = 1:length(MWDivKcats)
    %MWDivKcats_orig(i) = GetMWAndKcatForReaction(i, uniprots,kcats, swissprot, standardMW, swissProtHashMap);
    GetMWAndKcatForReaction(i, uniprots,kcats, swissprot, standardMW, swissProtHashMap);
    if rem(i,100) == 0 || i == length(MWDivKcats) 
        fprintf('.')
    end
end
toc%288.198142 with hash, the problem is reallocation of uniprots



compareVectorsNum(MWDivKcats, MWDivKcats_orig, '')
MWDivKcats2 = MWDivKcats;
MWDivKcats2(117) = NaN;
compareVectorsNum(MWDivKcats2, MWDivKcats_orig, '') %should fail, NaN, ok
MWDivKcats2(117) = 18;
compareVectorsNum(MWDivKcats2, MWDivKcats_orig, '') %should fail, numbers, ok

sel = find(MWDivKcats(1:100) ~= MWDivKcats_orig(1:100))
MWDivKcats(85)
MWDivKcats_orig(85)





%some test code (comparing to the ec model generated by Gecko):
%constructEquations(prepMinModel, 'HMR_3164No1') %gene rule: ENSG00000084754 and ENSG00000113790 and ENSG00000127884
%{'crotonyl-CoA[m] + H2O[m] + 4.9017e-08 prot_P30084[c] + 4.9017e-08 prot_P40939[c] + 4.9017e-08 prot_Q08426[c] => (S)-3-hydroxybutyryl-CoA[m]'}
%1/kcat = 4.9017e-08
%constructEquations(prepMinModel, 'draw_prot_P30084')
%{'31.3871 prot_pool[c] => prot_P30084[c]'}
%constructEquations(prepMinModel, 'draw_prot_P40939')
%{'82.9988 prot_pool[c] => prot_P40939[c]'}
%constructEquations(prepMinModel, 'draw_prot_Q08426')
%{'79.4942 prot_pool[c] => prot_Q08426[c]'}
% so, we would expect the outcome of MWDiv to be
%(31.3871+82.9988+79.4942)*4.9017e-08 %9.5034e-06
%ind = getIndexes(model, 'HMR_3164', 'rxns');

%MWDivKcat = GetMWAndKcatForReaction(model,ind, uniprots,kcats, swissprot, kegg, standardMW)%9.5034e-06, ok
%also check one with multiple enzymes
%HMR_4373 - gene rule ENSG00000105679 or ENSG00000111640
%constructEquations(prepMinModel, 'HMR_4373No1')
%{'pmet_HMR_4373[c] + 1.1871e-06 prot_O14556[c] => GAP[c] + NAD+[c] + Pi[c]'}
%constructEquations(prepMinModel, 'HMR_4373No2')
%{'pmet_HMR_4373[c] + 2.2474e-07 prot_P04406[c] => GAP[c] + NAD+[c] + Pi[c]'}
%constructEquations(prepMinModel, 'draw_prot_O14556')
%{'44.5008 prot_pool[c] => prot_O14556[c]'}
%constructEquations(prepMinModel, 'draw_prot_P04406')
%{'36.0528 prot_pool[c] => prot_P04406[c]'}
%alt 1: 
%44.5008*1.1871e-06%5.2827e-05
%alt 2:
%36.0528*2.2474e-07%8.1025e-06
%MWDivKcat2 = GetMWAndKcatForReaction(getIndexes(model, 'HMR_4373', 'rxns'), uniprots,kcats, swissprot, standardMW, swissProtHashMap)%8.1025e-06 - looks good
%constructEquations(model, 'HMR_4373') %looks good
%constructEquations(model, 'prot_pool_exchange')%looks good

%compare the resulting models before and after optimization (S matrix only)
cd C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy
load('data/ltModel.mat'); %this is expected to be an old one

ltModelNew = CreateECLightModel(ihuman);
ltModelNew = curateModel(ltModelNew);
ltModelNew.lb(ltModelNew.lb == -1000) = -Inf; %These operations help the solver, it runs faster with inf
ltModelNew.ub(ltModelNew.ub == 1000) = Inf;

bloodData = prepBloodData(false, true);
ltModelNewMin = minimizeModel(ltModelNew, bloodData);
ltModelNewPrep = setGrowthMedium(ltModelNewMin, true, 'Hams');


sel = ~contains(ltModel.rxns, ltModelNewPrep.rxns);
ltModel.rxns(sel)

compareVectorsNum(ltModelNewPrep.S, ltModel.S, 'S') %they usually differ on a few reactions, I don't think it matters much, probably some randomness there

%compareVectorsNum(ltModel.S, prepMinModelLite.S, 'S') %strange that these are slightly different (a few exchange reactions it seems)
%so, it seems that minimize model did something different. It does not
%matter much I'd say, so let's ignore it for now.


%look at the total pmets with only two reactions (draw + rxn), where we
%have filtered rev reactions that has a matching non-rev
%revRxnsSel = endsWith(prepMinModel.rxns, '_REV');
%revRxnsWithoutRev = extractBefore(prepMinModel.rxns(revRxnsSel),'_REV')
%toFilter = ismember(prepMinModel.rxns, revRxnsWithoutRev);
%sum(toFilter)
%STmp = prepMinModel.S(: ,~toFilter);
%numRxns2 = sum(STmp ~= 0, 2);
%numRxns2(protMets)
%toPlot = numRxns2(protMets);
%toPlot(toPlot > 20) = [];
%figure
%histogram(toPlot)


%numRxns(protMets)

