function [tumor,healthy] = ReadLung50000(pathchunk1, pathchunk2, pathchunk3, metaDataPath)
% ReadLung50000
%   Reads lung cancer data from file. Note that the file has been divided into 
%   chunks for matlab to be able to read it.
% Input:
%   pathchunk1          Path to the first data file
%   pathchunk2          Path to the second data file
%   pathchunk3          Path to the third data file
%   metaDataPath        Path to the metadata file
%
% Output: (optional)
%   tumor               SCDataset containing the cells from the tumors
%   healthy             SCDataset containing the cells from the healthy tissue
% Usage: [ds1,ds2] = ReadLung50000('../../ImportableData/LungCancer50000/AllCells_chunk1.txt', 
%                       '../../ImportableData/LungCancer50000/AllCells_chunk2.txt', 
%                       '../../ImportableData/LungCancer50000/AllCells_chunk3.txt', 
%                       '../../ImportableData/LungCancer50000/MetaData.txt');
%
% Johan Gustafsson, 2019-05-20
%

%path = 'C:/Work/MatlabCode/components/SCLib/ImportableData/LungCancer50000/TestData.txt';
%metaDataPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/LungCancer50000/MetaData.txt';

s = SCDataset;

%need to read it in chunks; requires too much memory otherwise
L = importdata(pathchunk1,'\t');
s.data = L.data;
s.cellIds = L.textdata(1, 2:end);
s.genes = L.textdata(2:end, 1);
L = importdata(pathchunk2,'\t');
s.cellIds = [s.cellIds L.textdata(1, 2:end)];
s.data = [s.data L.data];
L = importdata(pathchunk3,'\t');
s.cellIds = [s.cellIds L.textdata(1, 2:end)];
s.data = [s.data L.data];


%read and process metadata
t = readtable(metaDataPath, 'ReadVariableNames',true, 'ReadRowNames', false, 'Delimiter', '\t');
cellIds = table2cell(t(:,3));
cluster = table2cell(t(:,2));
cellType = table2cell(t(:,4));
patient = table2cell(t(:,9));
patient = cellfun(@(x) {num2str(x)}, patient);
fromTumor = strcmp(table2cell(t(:,6)),'TRUE');

[~,ia,ib] = intersect(s.cellIds, cellIds.');


numCells = size(cellIds,1);

temp = zeros(1,numCells);
%Get the cell type for all cells.
for i = 1:numCells
    clust = cluster(i);
    %for some subtypes another cell type should be set
    if strcmp(clust,'B_cell_4')
        temp(1,i) = Celltype.Mast;
    elseif strcmp(clust,'B_cell_8')
        temp(1,i) = Celltype.Dendritic;
    elseif strcmp(clust,'B_cell_9')
        temp(1,i) = Celltype.Erythroblast;
    elseif strcmp(clust,'Myeloid_5')
        temp(1,i) = Celltype.Langerhans;
    elseif strcmp(clust,'Myeloid_7')
        temp(1,i) = Celltype.Granulocyte;
    elseif strcmp(clust,'Myeloid_9') || strcmp(clust,'Myeloid_12')
        temp(1,i) = Celltype.Dendritic;
    elseif strcmp(clust,'T_cell_6')
        temp(1,i) = Celltype.NKCell;
    elseif strcmp(clust,'T_cell_2') || strcmp(clust,'T_cell_4') || strcmp(clust,'T_cell_5') || strcmp(clust,'T_cell_8') 
        temp(1,i) = Celltype.TCellCD8Pos;
    elseif strcmp(clust,'T_cell_1') || strcmp(clust,'T_cell_3') || strcmp(clust,'T_cell_9')
        temp(1,i) = Celltype.TCellCD4Pos;
    elseif strcmp(clust,'T_cell_7')
        temp(1,i) = Celltype.TCellReg;
    else
        %for the rest, just use the cell type specified
        ct = cellType(i);
        if strcmp(ct,'Alveolar')
            temp(1,i) = Celltype.Alveolar;
        elseif strcmp(ct,'B_cell')
            temp(1,i) = Celltype.BCell;
        elseif strcmp(ct,'EC')
            temp(1,i) = Celltype.Endothelial;
        elseif strcmp(ct,'Epi')
            temp(1,i) = Celltype.Epithelial;
        elseif strcmp(ct,'Fibro')
            temp(1,i) = Celltype.Fibroblast;
        elseif strcmp(ct,'Myeloid')
            temp(1,i) = Celltype.Macrophage;%took care of the other types above
        elseif strcmp(ct,'T_cell')
            temp(1,i) = Celltype.TCell;
        elseif strcmp(ct,'tumor')
            temp(1,i) = Celltype.Malignant;
        else
            temp(1,i) = Celltype.Unknown;
        end
    end
end
%now assign the cell types to the right rows, needs to be mapped by cell id!
s.cellType(1,ia.') = temp(ib).';

s = s.fillEmpties();%will create vectors of right size

%set the subclass if anyone wants to look at the clusters
temp = regexp(cluster, '.*_([0-9]+)', 'tokens');
strs = cellfun(@(c) c{1}{1},temp, 'UniformOutput', false);
nums = str2double(strs);
s.subCellType(1,ia.') = nums(ib).';

%set patient ids as sample ids - we lose some info about localization of
%sample within the tumor; can be set if desired
s.sampleIds(1,ia.') = patient(ib,1).';

%split in tumor and healthy
th(1,ia.') = fromTumor(ib,1).';
tumor = s.cellSubset(th);
healthy = s.cellSubset(~th);
%overwrite name with something more nice
tumor.name = 'lc 50000 tumor';
healthy.name = 'lc 50000 healthy';

end
