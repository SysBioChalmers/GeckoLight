classdef GeckoLightInstall
% GeckoInstall
%   Support for installing and uninstalling
%   Run GeckoLightInstall.install() to install (will set up the path in MATLAB)
%   Run GeckoLightInstall.uninstall() to clear the path from MATLAB
%
% Johan Gustafsson, 2021-06-11
%
    methods (Static)
        function install
            sourceDir = fileparts(which(mfilename));
            paths = GeckoLightInstall.GetFilteredSubPaths(sourceDir, '.*\.git.*');
            addpath(paths);
            savepath;
        end
        function uninstall
            sourceDir = fileparts(which(mfilename));
            paths = GeckoLightInstall.GetFilteredSubPaths(sourceDir, '.*\.git.*');
            rmpath(paths);
            savepath;
        end
        function path = getGeckoLightMainPath()
			path = fileparts(which(mfilename));
			path = strrep(path, '\', '/'); %get rid of backslashes in Windows
			if ~endsWith(path, '/')
				path = strcat(path,'/');
			end
		end

        function newPaths = GetFilteredSubPaths(path_, filter_)
            % Will fail if you have a directory containing ';'
            paths = genpath(path_);
            splitPaths = strsplit(paths, ';');
            %remove the last, it is empty
            splitPaths = splitPaths(1,1:end-1);
            matches = regexp(splitPaths, filter_, 'match');
            okPaths = cellfun(@isempty, matches);
            pathsLeft = splitPaths(1,okPaths);
            newPaths = strcat(char(join(pathsLeft,';')),';');
        end
    end
end
