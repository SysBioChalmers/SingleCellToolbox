classdef DsCFib
    % DsCFib
    %   Reads Cell line Fibroblast data from file into an SCDataset class
    %   Data is in counts. This is not a UMI dataset.
    %
    % Johan Gustafsson, 2020-01-23
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsCFib.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading CFibanoma...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/cfib.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsCFib.import('../../ImportableData/cfib/GSE66053_singlecellseq_star_HTSeqCounts.tsv');
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
        function  ds = import(path)
            % import
            %   Imports the data
            % Input:
            %   path       Path to the data file
            %
            % Usage: ds = DsCFib.import('../../ImportableData/cfib/GSE66053_singlecellseq_star_HTSeqCounts.tsv');
            %
            
            ds = SCDataset;
            ds.name = 'cfib';
            L = importdata(path,'\t');
            
            cellIdsTable = L.textdata(2:end,2);
            genesTable = L.textdata(2:end,3);
            dataTable = L.data(:,1);
            cellIds = unique(cellIdsTable);
            
            ds.cellIds = cellIds;
            sel = strcmp(cellIdsTable,cellIds(1));
            ds.genes = genesTable(sel);
            ds.data = zeros(length(ds.genes),length(ds.cellIds));
            for i = 1:length(cellIds)
                sel = strcmp(cellIdsTable,cellIds(i));
                ds.data(:,i) = dataTable(sel,:);
            end
            
            
            %genesToRem = {'__no_feature'; '__ambiguous'; '__too_low_aQual'; '__not_aligned'; '__alignment_not_unique'};
            %eh, just remove the 5 last genes...
            ds.data = ds.data(1:(end-5), :);
            ds.genes = ds.genes(1:(end-5), :);
            
            ds = ds.fillEmpties();
        end
    end
end

