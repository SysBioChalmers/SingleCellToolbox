classdef DsBC2
    % DsBC2
    %   Access wrapper for breast cancer data set BC2.
    %   The publication covers roughly 60000 cells, but only 47000 immune cells
    %   are included in this dataset; the rest need to be extracted from the raw files, which we have not done.
    %   Note that this code can read both raw and 'imputed' data, where the latter
    %   is treated with biscuit (see the article 'Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment')
    %   The imputed could potentially be better for an analysis where only this dataset
    %   is involved. We only read the raw data here
    %   The data are stored in files with the same format though.
    %
    % Johan Gustafsson, 2019-05-20
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsBC2.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading breast cancer 2...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/bc2.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsBC2.import('../../ImportableData/bc2_raw_corrected.csv', '../../ImportableData/bc2_cluster_ids.txt');
                    %convert to sparse matrix before saving not to exceed limit for
                    %saving. Also keep it that way.
                    v.data = sparse(v.data);
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
        function ds = import(path, clusterInfoPath)
            % import
            %   Import data from original files.
            % Input:
            %   path            Path to the data file,
            %   clusterInfoPath Path to the cell type info file
            %
            % Usage: ds = DsBC2.import('../../ImportableData/bc2_raw_corrected.csv', '../ImportableData/bc2_cluster_ids.txt');
            %
            
            %path = '../../ImportableData/bc2_raw_corrected.csv';
            %path = '../../ImportableData/bc2_test.csv';
            %clusterInfoPath = '../../ImportableData/bc2_cluster_ids.txt';
            
            %first read cluster info
            f = readtable(clusterInfoPath, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
            clusterTypesText = table2cell(f(:, 7));%creates a vertical cell array indexed by cluster id
            clusterCellTypes = cellfun(@DsBC2.String2CellTypeId, clusterTypesText,'UniformOutput',true);
            clusterDescription = table2cell(f(:, 6));
            
            ds = SCDataset;
            ds.name = 'scd_bc2';
            L = importdata(path,',');
            ds.data = L.data(:,4:end).';
            [m,n] = size(ds.data);
            ds.sampleIds = strcat(L.textdata(2:end, 1).','_',L.textdata(2:end, 2).');
            ds.cellIds = arrayfun(@(x) int2str(x), L.data(:,3),'UniformOutput',false);
            ds.cellIds = ds.cellIds.';
            ds.genes = L.textdata(1, 6:end).';
            ds.cellType = arrayfun(@(x) clusterCellTypes(x),L.data(:,2));
            ds.cellType = ds.cellType.';
            ds.extraCellInfo = arrayfun(@(x) clusterDescription(x),L.data(:,2));
            ds.extraCellInfo = ds.extraCellInfo.';
            ds = ds.fillEmpties();
        end
        
        function ret = String2CellTypeId(str)
            if strcmp(str,'Tcell')
                ret = Celltype.TCell;
            elseif strcmp(str,'Tcell_CD4')
                ret = Celltype.TCellCD4Pos;
            elseif strcmp(str,'Tcell_CD8')
                ret = Celltype.TCellCD8Pos;
            elseif strcmp(str,'Tcell_Reg')
                ret = Celltype.TCellReg;
            elseif strcmp(str,'Bcell')
                ret = Celltype.BCell;
            elseif strcmp(str,'NKcell')
                ret = Celltype.NKCell;
            elseif strcmp(str,'DC')
                ret = Celltype.Dendritic;
            elseif strcmp(str,'Macrophage')
                ret = Celltype.Macrophage;
            elseif strcmp(str,'Monocyte')
                ret = Celltype.Monocyte;
            elseif strcmp(str,'Neutrophil')
                ret = Celltype.Neutrophil;
            elseif strcmp(str,'Mast')
                ret = Celltype.Mast;
            elseif strcmp(str,'Unknown')
                ret = Celltype.Unknown;
            else
                error(strcat('Unknown type: ', str));
                ret = Celltype.Unknown;
            end
        end
    end
end
