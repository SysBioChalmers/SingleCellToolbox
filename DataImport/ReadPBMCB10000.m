function ds = ReadPBMCB10000(directoryPath)
% ReadPBMCB10000
%   Reads the PBMCB10k cells from file into an SCDataset.  
%   The dataset covers roughly 10000 B cells from blood.
% Input:
%   directoryPath       Path to the 10x files. No slash at the end.
%
% Usage: ds = ReadPBMCB10000('../../ImportableData/PBMC10000BCells/filtered_matrices_mex/hg19');
%
% Johan Gustafsson, 2019-05-20
%

%directoryPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/PBMC10000BCells/filtered_matrices_mex/hg19';

ds = Read10xMatrix(directoryPath);
ds.cellType(1,:) = Celltype.BCell;
ds.name = 'pbmc b 10000';

end
