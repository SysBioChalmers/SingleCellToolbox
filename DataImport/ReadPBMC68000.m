function ds = ReadPBMC68000(directoryPath)
% ReadPBMC68000
%   Reads the PBMC68k cells from file into an SCDataset.  
%   The dataset covers roughly 68000 cells from blood.
%   The file 68k_pbmc_barcodes_annotation.tsv is expected in directoryPath
% Input:
%   directoryPath       Path to the 10x files. No slash at the end.
%
% Usage: ds = ReadPBMC68000('../../ImportableData/PBMC68000PatAFresh/filtered_matrices_mex/hg19');
%
% Johan Gustafsson, 2019-05-20
%


%directoryPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/PBMC68000PatAFresh/filtered_matrices_mex/hg19';
classificationPath = strcat(directoryPath,'/68k_pbmc_barcodes_annotation.tsv');

ds = Read10xMatrix(directoryPath);
ds.name = 'pbmc 68000';

f = readtable(classificationPath, 'ReadRowNames', false, 'Delimiter', '\t', 'FileType', 'text');
c = table2cell(f);
ct = cellfun(@String2CellTypeId, c(:,4));

[~, ia, ib] = intersect(ds.cellIds.', c(:, 3));
ds.cellType(1,ia) = ct(ib,1).';


    function ret = String2CellTypeId(str)
        %we ignore all detailed subsets for now
        if strcmp(str,'CD8+ Cytotoxic T')
            ret = Celltype.TCellCD8Pos;
        elseif strcmp(str,'CD8+/CD45RA+ Naive Cytotoxic')
            ret = Celltype.TCellCD8Pos;
        elseif strcmp(str,'CD4+/CD45RO+ Memory')
            ret = Celltype.TCellCD4Pos;
        elseif strcmp(str,'CD4+/CD45RA+/CD25- Naive T')
            ret = Celltype.TCellCD4Pos;
        elseif strcmp(str,'CD4+ T Helper2')
            ret = Celltype.TCellCD4Pos;
        elseif strcmp(str,'CD4+/CD25 T Reg')
            ret = Celltype.TCellReg;
        elseif strcmp(str,'CD19+ B')
            ret = Celltype.BCell;
        elseif strcmp(str,'CD14+ Monocyte')
            ret = Celltype.Monocyte;
        elseif strcmp(str,'CD56+ NK')
            ret = Celltype.NKCell;
        elseif strcmp(str,'Dendritic')
            ret = Celltype.Dendritic;
        elseif strcmp(str,'CD34+')
            ret = Celltype.HematopeticStemOrProgenitor;
        %elseif strcmp(str,'Megakaryocytes')
        %    ret = Celltype.Megakaryocyte;
        else
            disp(strcat('unknown type: ',str));
            ret = Celltype.Unknown;
        end
    end

end
