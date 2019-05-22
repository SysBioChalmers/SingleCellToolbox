function ds = ReadGSE112845(path, classificationPath)
% ReadGSE112845
%   Reads one of the GSE112845 datasets from file into an SCDataset. Note 
%   that there are several datasets to choose from.
% Input:
%   path                Path to the 10x data folder
%   classificationPath  Path to the cell type info file, can be omitted
%
% Usage: ds = ReadGSE112845('../../ImportableData/GSE112845/DTM-X_PBMC_live', '../../ImportableData/GSE112845/DTM-X_PBMC_live_ct.txt');
%
% Johan Gustafsson, 2019-05-20
%

%path = 'C:/Work/MatlabCode/components/SCLib/ImportableData/GSE112845/DTM-X_PBMC_live';
%classificationPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/GSE112845/DTM-X_PBMC_live_ct.txt';

ds = Read10xMatrix(path);
ds.name = 'GSE112845';

if ~isempty(classificationPath)
    %read ds.cellType file from the authors
    f = readtable(classificationPath, 'ReadRowNames', false, 'Delimiter', '\t');
    c = table2cell(f);
    %transform from string type to SCLib cell type enum
    ct = cellfun(@String2CellTypeId, c(:,2));

    [~, ia, ib] = intersect(ds.cellIds, c(:, 1));
    ds.cellType(1,ia) = ct(ib,1).';

    %remove the cells with cell type 'unknown' - these are low quality cells
    %that did not have a classification in the file
    ds = ds.cellSubset(ds.cellType ~= Celltype.Unknown);
else
    %assume it is the CD8+ dataset
    ds.cellType(1,:) = Celltype.TCellCD8Pos;
    ds.cellType = ds.cellType;
end


    function ret = String2CellTypeId(str)
        if strcmp(str,'CD4 T cells')
            ret = Celltype.TCell;%subclass of t cells not reliable
        elseif strcmp(str,'CD8 T cells')
            ret = Celltype.TCell;
        elseif strcmp(str,'B cells')
            ret = Celltype.BCell;
        elseif strcmp(str,'CD14+ Monocytes')
            ret = Celltype.Monocyte;%ignore subset for now
        elseif strcmp(str,'FCGR3A+ Monocytes')
            ret = Celltype.Monocyte;
        elseif strcmp(str,'NK cells')
            ret = Celltype.NKCell;
        elseif strcmp(str,'Dendritic cells')
            ret = Celltype.Dendritic;
        elseif strcmp(str,'Megakaryocytes')
            ret = Celltype.Megakaryocyte;
        else
            disp(strcat('unknown type: ',str));
            ret = Celltype.Unknown;
        end
    end
end
