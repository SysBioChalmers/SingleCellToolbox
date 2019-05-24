classdef SingleCellToolboxInstall
% DSAVEInstall
%   Support for installing and uninstalling DSAVE
%   Run DSAVEInstall.install() to install (will set up the path in MATLAB)
%   Run DSAVEInstall.uninstall() to clear the path from MATLAB
%
% Johan Gustafsson, 2019-05-20
%
    methods (Static)
        function install
            sourceDir = fileparts(which(mfilename));
            paths = SingleCellToolboxInstall.GetFilteredSubPaths(sourceDir, '.*\.git.*');
            addpath(paths);
            savepath;
        end
        function uninstall
            sourceDir = fileparts(which(mfilename));
            paths = SingleCellToolboxInstall.GetFilteredSubPaths(sourceDir, '.*\.git.*');
            rmpath(paths);
            savepath;
        end
        
        function newPaths = GetFilteredSubPaths(path_, filter_)
            % Will fail if you have a directory containing ';'
            paths = genpath(path_);
            splitPaths = strsplit(paths, ';');
            %remove the last, it is empty
            splitPaths = splitPaths(1,1:end-1);
            matches = regexp(splitPaths, '.*\.git.*', 'match');
            okPaths = cellfun(@isempty, matches);
            pathsLeft = splitPaths(1,okPaths);
            newPaths = strcat(char(join(pathsLeft,';')),';');
        end
        
    end
end
