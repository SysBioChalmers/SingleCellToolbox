classdef DsB10k
    % DsB10k
    %   Reads the DsB10k cells from file into an SCDataset.
    %   The dataset covers roughly 10000 B cells from blood.
    %
    % Johan Gustafsson, 2019-05-22
    %
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsB10k.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading pbmcb10000 data ...');
                filename = '../../TempData/pbmcb10000.mat';
                prevDir = DsHelper.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsB10k.import('../../ImportableData/PBMC10000BCells/filtered_matrices_mex/hg19');
                    save(filename, 'v');
                else
                    a = load(filename);
                    v = a.v;
                end
                DsHelper.restoreDir(prevDir);
            end
            ret = v;
        end
    end
    
    methods(Static, Access = private)
        function ds = import(directoryPath)
            % import
            %   Imports the data
            % Input:
            %   directoryPath       Path to the 10x files. No slash at the end.
            %
            % Usage: ds = DsB10k.import('../../ImportableData/PBMC10000BCells/filtered_matrices_mex/hg19');
            %
            
            %directoryPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/PBMC10000BCells/filtered_matrices_mex/hg19';
            
            ds = Read10xMatrix(directoryPath);
            ds.cellType(1,:) = Celltype.BCell;
            ds.name = 'pbmc b 10000';
            
        end
    end
end