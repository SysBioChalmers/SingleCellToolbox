function ds = ReadPBMCTCD4Mem10000(directoryPath)
% ReadPBMCTCD4Mem10000
%   Reads the PBMCTCDMem10k cells from file into an SCDataset.  
%   The dataset covers roughly 10000 CD4+ memory T cells from blood.
%   This dataset comes from the publication: Zheng et al, 
%   “Massively parallel digital transcriptional profiling of single cells”
% Input:
%   directoryPath       Path to the 10x files. No slash at the end.
%
% Usage: ds = ReadPBMCTCD4Mem10000('../../ImportableData/PBMCCD4TCellsMemory/filtered_matrices_mex/hg19');
%
% Johan Gustafsson, 2019-05-20
%

ds = Read10xMatrix(directoryPath);
ds.cellType(1,:) = Celltype.TCellCD4Memory;
ds.name = 'pbmc T CD4Mem 10000';

end
