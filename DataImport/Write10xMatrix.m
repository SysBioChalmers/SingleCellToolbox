%for testing:
%directoryPath = 'data/exportToR/datasets/TCD8'
%[~,~,ds] = DsGSE112845.get();
function ds = Write10xMatrix(directoryPath, ds)
% Write10xMatrix
%   Writes a SCDataset to the 10x Genomics Matrixmarket format
%
% Input:
%   directoryPath    Path to the directory in which the 10x files are,
%                    i.e. the files  'matrix.mtx', 'genes.tsv' and 
%                    'barcodes.tsv'. No slash at the end!
%   ds               Dataset
%
% Usage: ds = Write10xMatrix('../../ImportableData/SomePath', ds)
%
% Johan Gustafsson, 2021-10-14
%

if ~exist(directoryPath, 'dir')
   mkdir(directoryPath)
end

barcodesFilename = 'barcodes.tsv';
genesFilename = 'genes.tsv';
matrixFilename = 'matrix.mtx';
matrixPath = strcat(directoryPath, '/', matrixFilename);
genesPath = strcat(directoryPath, '/', genesFilename);
barcodesPath = strcat(directoryPath, '/', barcodesFilename);

%write genes
fileID = fopen(genesPath, 'w');
for i = 1:length(ds.genes)
    fprintf(fileID, '%s\n', ds.genes{i});
end
fclose(fileID);

%write barcodes
fileID = fopen(barcodesPath, 'w');
for i = 1:length(ds.cellIds)
    fprintf(fileID, '%s\n', ds.cellIds{i});
end
fclose(fileID);

%now write the matrix (probably slow)
%header:
fileID = fopen(matrixPath, 'w');
fprintf(fileID, '%%%%MatrixMarket matrix coordinate real general\n%%\n');
%matrix size row:
fprintf(fileID, '%d %d %d\n', length(ds.genes), length(ds.cellIds), full(sum(sum(ds.data ~= 0, 1),2)));
%the matrix
[i,j,s] = find(ds.data);
fprintf(fileID, '%d %d %d\n', [i j s].');
%finished
fclose(fileID);

end
