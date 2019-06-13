%Reads 10x files into a SCDataset, 
function ds = Read10xMatrix(directoryPath, matrixFilename, genesFilename, ...
                            barcodesFilename, barcodesHeaders, genesColumn)
% Read10xMatrix
%   Reads the standard output single-cell data from 10x Genomics into an
%   SCDataset
%
% Input:
%   directoryPath    Path to the directory in which the 10x files are,
%                    i.e. the files  'matrix.mtx', 'genes.tsv' and 
%                    'barcodes.tsv'. No slash at the end!
%                    The barcodes file may contain multiple columns as long
%                    as the barcodes are in the first!
%   matrixFilename   (Optional) Matrix filename, defaults to 'matrix.mtx'
%   genesFilename    (Optional) Genes filename, defaults to 'matrix.mtx'
%   barcodesFilename (Optional) Barcodes filename, defaults to 'matrix.mtx'
%   barcodesHeaders  (Optional) If the barcodes file contains headers,
%                    defaults to false
%   genesColumn      (Optional) Column for the genes in the genes file,
%                    defaults to 2
%
% Usage: ds = Read10xMatrix('../../ImportableData/SomePath')
%
% Johan Gustafsson, 2019-06-02
%

if nargin < 6
     genesColumn = 2;
end
if nargin < 5
     barcodesHeaders = false;
end
if nargin < 4
     barcodesFilename = 'barcodes.tsv';
end
if nargin < 3
     genesFilename = 'genes.tsv';
end
if nargin < 2
     matrixFilename = 'matrix.mtx';
end



matrixPath = strcat(directoryPath, '/', matrixFilename);
genesPath = strcat(directoryPath, '/', genesFilename);
barcodesPath = strcat(directoryPath, '/', barcodesFilename);

ds = SCDataset;

%couldn't get readmtx to work, using fscanf in combination with spconvert
%instead
%first count the number of comment rows
fileID = fopen(matrixPath);
numCommentRows = 0;
while 1
    tline = fgetl(fileID);
    if (tline(1) == '%')
        numCommentRows = numCommentRows + 1;
    else
        break;
    end
end
fclose(fileID);

%now start reading the file again
fileID = fopen(matrixPath);
%get rid of the comment rows
for i = 1:numCommentRows
    fgetl(fileID);
end

%get rid of size
fgetl(fileID);
%now read
formatSpec = '%d\t%d\t%d';
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
A = A.';
ds.data = spconvert(A);

fclose(fileID);

f = readtable(genesPath, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t', 'FileType', 'text');
ds.genes = table2cell(f(:, genesColumn));

%fill out the last few genes with zeros, they are not part of the matrix
[matrows,cols] = size(ds.data);
generows = size(ds.genes,1);
if (matrows < generows)
    ds.data = [ds.data;zeros(generows-matrows,cols)];
end

f = readtable(barcodesPath, 'ReadVariableNames', barcodesHeaders, 'ReadRowNames', false, 'Delimiter', '\t', 'FileType', 'text');
ds.cellIds = table2cell(f(:, 1)).';

ds = ds.fillEmpties();

end
