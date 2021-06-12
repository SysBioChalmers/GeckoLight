%Abstract Base class for adapters for different species
classdef (Abstract) SpeciesAdapter 
	methods (Abstract)
		parameters = getParameters(obj)
        result = getFilePath(obj, filename)
		[spont,spontRxnNames] = getSpontaneousReactions(obj,model);
    end
	methods
		function model = manualModifications(obj,model) %default is to do nothing
		end
	end
end
