%Reads breast cancer data from file



function ds = ReadBC2(path, clusterInfoPath)
% ReadBC2
%   Reads breast cancer data from file. 
%   The publication covers roughly 60000 cells, but only 47000 immune cells
%   are included in this dataset; the rest need to be extracted from the raw files.
%   Note that this code can read both raw and 'imputed' data, where the latter
%   is treated with biscuit (see the article 'Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment')
%   The imputed one is probably better for an analysis where only this dataset
%   is involved.
%   The data are stored in files with the same format.%
% Input:
%   path            Path to the data file,
%   clusterInfoPath Path to the cell type info file
%
% Usage: ds = ReadBC2('../../ImportableData/bc2_raw_corrected.csv', '../ImportableData/bc2_cluster_ids.txt');
%
% Johan Gustafsson, 2019-05-20
%




%path = '../../ImportableData/bc2_raw_corrected.csv';
%path = '../../ImportableData/bc2_test.csv';
%clusterInfoPath = '../../ImportableData/bc2_cluster_ids.txt';

%first read cluster info
f = readtable(clusterInfoPath, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
clusterTypesText = table2cell(f(:, 7));%creates a vertical cell array indexed by cluster id
clusterCellTypes = cellfun(@String2CellTypeId, clusterTypesText,'UniformOutput',true);
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
