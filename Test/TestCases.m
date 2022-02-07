load('C:/Work/MatlabCode/components/human-GEM/Human-GEM/model/Human-GEM.mat')
cd 'C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy'
speciesAdapter = HumanGEMAdapter();
minAcceptableKCat = 0;

parameters = speciesAdapter.getParameters();

oldMetNames = model.metNames;

%Remove blocked rxns + correct model.rev:
[model,name,modelVer] = preprocessModelLt(ihuman,'ecLight','1');


%%%%%%%%%%%%%%%%%%
%GetEnzymeCodes
%%%%%%%%%%%%%%%%%%


%cd C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy/GeckoWork
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
[KCATcell, SAcell] = loadBRENDAdataLt(speciesAdapter);
%[KCATcell, SAcell] = loadBRENDAdata;
%cd C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy/GeckoWork
tic
%kcats      = matchKcatsOpt(model_data,parameters.org_name, KCATcell, SAcell, 'C:/Work/MatlabCode/components/Gecko/GeckoForHuman20210603/GECKO/Databases/PhylDist.mat');
kcats = matchKcatsOpt(model_data,parameters.org_name, KCATcell, SAcell, speciesAdapter, minAcceptableKCat);
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

%Some code below is copied from CreateECLightModel

%Get kcat value for both directions:
Fkcat = kcats.forw.kcats;
Bkcat = kcats.back.kcats;
rev   = logical(model_data.model.rev);
kcatVals = [Fkcat;Bkcat(rev,:)];
wcLevels = [kcats.forw.wcLevel;kcats.back.wcLevel(rev,:)];

%Update uniprots with both directions:
uniprots = [model_data.uniprots; model_data.uniprots(rev,:)];

%Update matched genes with both directions:
matchedGenes = [model_data.matchedGenes; model_data.matchedGenes(rev,:)];

%Convert to irreversible model with RAVEN function (will split in 2 any reversible rxn):
irrevModel = convertToIrrev(model_data.model);


%What we want to find for each reaction
%is the min value of MW/kcat when comparing all parallel reactions (i.e. OR in the GR rules)
%The plan is then to add one reaction with prot_pool_exchange, one
%metabolite called prot_pool, and a stochiometric coefficient of MW/kcat of
%that for each reaction. For complexes, the MW are summed up, since the kcat is the same.

%Load databases:
data      = load(speciesAdapter.getFilePath('ProtDatabase.mat')); %this is loaded twice, (also in a function above), a bit unnecessary
swissprot = data.swissprot;

standardMW = median([swissprot{:,5}])/1000;%the median of all the proteins

%Compare a few reactions to see that they have been assigned the right values.
%We assume that the standard Gecko works, which means that we can trust the values of the MW and kcat above
%look at a few reactions, one with a simple GPR, one with AND, and one with OR

%First test: simple GPR
%MAR04097: acetate[c] + ATP[c] + CoA[c] => acetyl-CoA[c] + AMP[c] + PPi[c]  ENSG00000131069

%Generate a Gecko Light model:
glModel = CreateECLightModel(ihuman, true, 0, speciesAdapter);


ind1 = find(strcmp(irrevModel.rxns,'MAR04097'));
kcat1 = kcatVals(ind1,1) %there is only one, 21600
protNames1 = uniprots{ind1,1};
MW1 = swissprot{strcmp(protNames1,swissprot(:,1)),5}/1000;
prodVal1 = -glModel.S(length(glModel.mets), ind1);
gt1 = MW1/kcat1;

%check:
gt1 == prodVal1 %ok

%second test: AND rule
%MAR03208: ATP[m] + HCO3-[m] + propanoyl-CoA[m] => ADP[m] + H+[m] + methylmalonyl-CoA[m] + Pi[m]     ENSG00000114054 and ENSG00000175198

ind2 = find(strcmp(irrevModel.rxns,'MAR03208'));
kcat2 = kcatVals(ind2,1); %there is only one, it is a complex
protNames2 = uniprots{ind2,1:10};
protNames2spl = strsplit(protNames2);
MW2_1 = swissprot{strcmp(protNames2spl{1},swissprot(:,1)),5}/1000;
MW2_2 = swissprot{strcmp(protNames2spl{2},swissprot(:,1)),5}/1000;
prodVal2 = full(-glModel.S(length(glModel.mets), ind2));
gt2 = (MW2_1+MW2_2)/kcat2;

%check:
abs(gt2 - prodVal2) < 10^-10 %ok


%third test: OR rule
%MAR01987:17alpha-hydroxyprogesterone[r] + NADPH[r] + O2[r] <=> 4-androstene-3,17-dione[r] + acetate[r] + H2O[r] + NADP+[r]     ENSG00000148795 or ENSG00000197580

ind3 = find(strcmp(irrevModel.rxns,'MAR01987'));
kcats3 = kcatVals(ind3,1:2);
kcats3
protNames3 = uniprots(ind3,1:2);
MW3_1 = swissprot{strcmp(protNames3{1},swissprot(:,1)),5}/1000;
MW3_2 = swissprot{strcmp(protNames3{2},swissprot(:,1)),5}/1000;
prodVal3 = full(-glModel.S(length(glModel.mets), ind3));
gt3_1 = MW3_1/kcats3(1);
gt3_2 = MW3_2/kcats3(2);
gt3 = min(gt3_1,gt3_2);

%check:
abs(gt3 - prodVal3) < 10^-10 %ok



%Some old test code below, not used
%%%%%%%%%%%%%%%%%%%%%%%%%


%tic
%[MWDivKcats,rxnWcLevels] = GetMWAndKcats(uniprots,kcatVals,wcLevels, swissprot, standardMW);
%toc% Elapsed time is 2.533362 seconds.

%compare


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
%cd C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Johan/OptimalTMEGrowthStrategy
%load('data/ltModel.mat');
%load('data/prepMinModelLite.mat');
%cd GeckoWork

%sel = ~contains(ltModel.rxns, prepMinModelLite.rxns);
%ltModel.rxns(sel)

%compareVectorsNum(ltModelOrig.S, ecModelOrigLite.S, 'S') %ok
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

