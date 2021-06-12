%tc0001 - getSpontaneous
%{
load(strcat(GetHumanGEMRootPath(),'model/Human-GEM.mat'));

tic
rxns_tsv = importTsvFile(strcat(GetHumanGEMRootPath(),'model/reactions.tsv'));
spontTmp = rxns_tsv.spontaneous;
spontRxnNames = rxns_tsv.rxns;
spontaneous1 = NaN(length(ihuman.rxns),1);
%simple but slow
for i = 1:length(ihuman.rxns)
   spontaneous1(i) = spontTmp(strcmp(spontRxnNames, ihuman.rxns(i)));
end
toc %Elapsed time is 8.435539 seconds.

spAd = HumanGEMAdapter();
tic
spontaneous2 = spAd.getSpontaneousReactions(ihuman);
toc %Elapsed time is 1.910023 seconds. - includes loading
all(spontaneous1 == spontaneous2) %ok
%}