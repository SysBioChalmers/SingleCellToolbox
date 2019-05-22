classdef DsTCD4Mem
    % DsTCD4Mem
    %   Reads the DsTCD4Mem cells from file into an SCDataset.
    %   The dataset covers roughly 10000 CD4+ memory T cells from blood.
    %   This dataset comes from the publication: Zheng et al,
    %   “Massively parallel digital transcriptional profiling of single cells”
    %
    % Johan Gustafsson, 2019-05-22
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsTCD4Mem.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading pbmccd4mem10000 data ...');
                filename = '../../TempData/pbmctcd4mem10000.mat';
                prevDir = DsHelper.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsTCD4Mem.import('../../ImportableData/PBMCCD4TCellsMemory/filtered_matrices_mex/hg19');
                    save(filename, 'v');
                else
                    a = load(filename);
                    v = a.v;
                end
                DsHelper.restoreDir(prevDir);
            end
            ret = v;
        end
        
        
        
        function ds = import(directoryPath)
            % import
            %   Imports the data
            % Input:
            %   directoryPath       Path to the 10x files. No slash at the end.
            %
            % Usage: ds = DsTCD4Mem.import('../../ImportableData/PBMCCD4TCellsMemory/filtered_matrices_mex/hg19');
            %
            
            ds = Read10xMatrix(directoryPath);
            ds.cellType(1,:) = Celltype.TCellCD4Memory;
            ds.name = 'pbmc T CD4Mem 10000';
            
        end
    end
end