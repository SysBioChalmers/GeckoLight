function parameters = getModelParameters
% getModelParameters
%
%   Set model and organism specific parameters that are used by the
%   ecModel generation pipeline.
%
%   Hao Wang Last edited: 2020-08-23
%

%Average enzyme saturation factor
parameters.sigma = 0.5;

%Total protein content in the cell [g protein/gDw]
parameters.Ptot = 0.505717472;  %Average across NCI60 cell lines

%Minimum growth rate the model should grow at [1/h]
parameters.gR_exp = 1E-3; %[g/gDw h]/Average across NCI60 cell lines

%Provide your organism scientific name
parameters.org_name = 'rattus norvegicus';

%Provide your organism KEGG ID
parameters.keggID = 'rno';

%The name of the exchange reaction that supplies the model with carbon (rxnNames)
parameters.c_source = 'MAR09034'; 

%Rxn Id for biomass pseudoreaction
parameters.bioRxn = 'MAR00021';

%Experimental carbon source uptake (optional)
parameters.c_UptakeExp = 1; %[mmol/gDw h]/Average across NCI60 cell lines

%Compartment name in which the added enzymes should be located
parameters.enzyme_comp = 'Cytosol';
end