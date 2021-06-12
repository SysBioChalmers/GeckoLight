classdef HumanGEMAdapter < SpeciesAdapter 
    methods
		function parameters = getParameters(obj)
			%these parameters are just copied from getParams in ecModels
			%Average enzyme saturation factor
			parameters.sigma = 0.5;

			%Total protein content in the cell [g protein/gDw]
			parameters.Ptot = 0.505717472;  %Average across NCI60 cell lines

			%Minimum growth rate the model should grow at [1/h]
			parameters.gR_exp = 0.020663429; %[g/gDw h]/Average across NCI60 cell lines

			%Provide your organism scientific name
			parameters.org_name = 'homo sapiens';

			%Provide your organism KEGG ID
			parameters.keggID = 'hsa';

			%The name of the exchange reaction that supplies the model with carbon (rxnNames)
			parameters.c_source = 'HMR_9034'; 

			%Experimental carbon source uptake (optional)
			parameters.c_UptakeExp = 0.641339301; %[mmol/gDw h]/Average across NCI60 cell lines

			%Rxn Id for non-growth associated maitenance pseudoreaction
			parameters.NGAM = 'biomass_maintenance_Recon3D'; %this is not used

			%Compartment name in which the added enzymes should be located
			parameters.enzyme_comp = 'Cytosol';

		end
		
        function result = getFilePath(obj, filename)
			result = strcat(GeckoLightInstall.getGeckoLightMainPath(), 'data/humanGEM/', filename);
		end
		
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			rxns_tsv = importTsvFile(strcat(GetHumanGEMRootPath(),'model/reactions.tsv'));
			spont = rxns_tsv.spontaneous;
			spontRxnNames = rxns_tsv.rxns;
		end
		
		function model = manualModifications(obj,model) %default is to do nothing
			%So, there are some of these in ecModels - it is a bit unclear if any of these are relevant here
			%we do nothing for now.
			%In general, manual modifications should be done to the model before sending it in.
		end
	end
end
