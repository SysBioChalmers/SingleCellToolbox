function ds = DSAVEAlignDataset(inDs, templInfo, progrBarCtxt)
% DSAVEAlignDataset
%   Aligns the dataset to the template, by removing genes, removing cells and
%   down-sampling the data to match the template as good as possible.
%   All datasets successfully aligned to the same template will have almost
%   identical sampling noise.
% Input:
%   inDs            The input dataset (cell population) to align.
%   templInfo       The template to align to
%   progrBarCtxt    (optional) Progress bar context.
%
% Usage: ds = DSAVEAlignDataset(inDs, templInfo)
%
% Johan Gustafsson, 2019-05-20
%

if nargin < 3
    progrBarCtxt = ProgrBarContext;
end

progbar = ProgrBar(['Aligning dataset to template ''' inDs.name ''''], progrBarCtxt);

%create the adapted dataset of the right size and with the right genes
inDs = inDs.randSample(size(templInfo.UMIDistr, 2));
inDs = inDs.geneSubset(templInfo.geneSet);
ds = inDs;
ds = ds.geneSubset(templInfo.geneSet);
ds.name = ['Aligned dataset from ds ' inDs.name];
numGenes = size(ds.genes,1);
numCells = size(ds.data,2);
%zero the data
ds.data = zeros(numGenes,numCells);%sparsify in the end

%then do downsampling
origUMIs = full(sum(inDs.data,1));%sum of UMIs for each cell
matchUMIs = templInfo.UMIDistr;%sum of UMIs for each cell
%we need create matching UMIs from the template. The template does not
%necessarily have the same number of cells, so we need to handle that.
%First duplicate the template until there are fewer cells left than in the
%template; then randomly select from the template to fill out the rest.

%{ 
%This code was for if we wanted to use all cells in a dataset. We concluded
in the end that this will bias the results and that it is not good to do.
ts = size(templInfo.UMIDistr,2);
numDupl = floor(numCells/ts);
cellsToRand = numCells - ts*numDupl;

matchUMIs = zeros(1,numCells);%sum of UMIs for each cell
%duplicate (if enough cells)
for i = 1:numDupl
    matchUMIs(1,(1+(i-1)*ts):i*ts) = templInfo.UMIDistr;
end
%fill randomly in the end
if (cellsToRand ~= 0)
    matchUMIs(1,(1+numDupl*ts):end) = templInfo.UMIDistr(1,randsample(ts,cellsToRand));
end
%}
%sort both UMI vectors and try to downsample to match the match set
%as good as possible
[~,iOrig] = sort(origUMIs);
[~,iMatchUMIs] = sort(matchUMIs);
idealUMIs = zeros(1,numCells);
idealUMIs(iOrig) = matchUMIs(iMatchUMIs);
tooSmall = origUMIs < idealUMIs;
newUMIs = min(origUMIs, idealUMIs);
UMIsLost = sum(idealUMIs(tooSmall)- newUMIs(tooSmall));
UMIsToSpend = UMIsLost;
toRemUMIs = origUMIs - newUMIs;
%for some cells, the UMIs 
%spread the lost UMIs over the other cells
i = 1;
%take care of the fact that the match dataset may
%have more UMIs than the ds to downsample
if (UMIsToSpend > sum(toRemUMIs))
    disp(strcat('Warning: Failed to downsample due to that the template had more UMIs than the target dataset!'));
    toRemUMIs = zeros(1,numCells);
else

    while UMIsToSpend > 0 
        if (toRemUMIs(1,i) > 0)
            toRemUMIs(1,i) = toRemUMIs(1,i) - 1;
            UMIsToSpend = UMIsToSpend - 1;
        end
        i = i+1;
        if i > numCells
            i = 1;
        end
    end
end

%Loop through all cells
for i = 1:numCells
    progbar.Progress(i/numCells);
    %first, randomly select a number of indexes from the total UMIs to
    %remove
    indexesToRem = randperm(origUMIs(1,i), toRemUMIs(1,i));
    %Then create index edges for each gene, so if a gene has 5 UMIs the
    %edges to the left and right will differ 5.
    cellData = inDs.data(:,i);
    edges = [0;cumsum(cellData)+0.1];%we need to add 0.1 to the edges since histcounts checks if edge(k) <= index < edge(k+1). Otherwise, index 1 would not end up between 0 and 1
    %Now get the number of index hits you get within the edge range for
    %each gene
    [subtr,~] = histcounts(indexesToRem,edges);
    %size(inDs.data(:,i))
    %size(subtr)
    ds.data(:,i) = inDs.data(:,i) - subtr.';
end

%now sparsify, will be slower to work with a sparse matrix in the loop
ds.data = sparse(ds.data);

progbar.Done();

end

