function  ds = ReadLiverTCells(pathData, pathCellFilter)
% ReadLiverTCells
%   Reads Liver cancer T cells from file into an SCDataset.  
%   The dataset covers roughly 4000 cells.
%   No classification other than that they are T cells is available.
% Input:
%   path                Path to the data file
%   pathCellFilter      Path to the cell quality file
%
% Usage: ds = ReadLiverTCells('../../ImportableData/GSE98638_HCC.TCell.S5063.count.txt', '../../ImportableData/GSE98638_HCC.TCell.OkCellIds.txt');
%
% Johan Gustafsson, 2019-05-20
%

ds = SCDataset;
ds.name = 'LCTCells';
L = importdata(pathData,'\t');
ds.data = L.data;
[m,n] = size(ds.data);
ds.cellIds = L.textdata(1, 3:end);
ds.genes = L.textdata(2:end, 2);
%filter on cell filter, i.e. remove the cells that did not pass quality
%control in the paper:
f = readtable(pathCellFilter, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
okCells = table2cell(f);
[ds.cellIds, ia, ib] = intersect(ds.cellIds, okCells);
ds.data = ds.data(:, ia);

%extract patient id, format is either 'xxxx-xx-yyyy' or 'xxxx-yyyy', and we
%only want the y:s
%a bit annoying, we have to extract two tokens and only keep the last
temp = regexp(ds.cellIds, '\w+-(\w+-)?(\w+)', 'tokens');
ds.sampleIds = cellfun(@(c) c{1}{2},temp, 'UniformOutput', false);
ds.cellType = repmat(Celltype.TCell,size(ds.genes,1),size(ds.cellIds,2));

ds = ds.fillEmpties();
