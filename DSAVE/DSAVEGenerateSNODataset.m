function ds = DSAVEGenerateSNODataset(templDs, progrBarCtxt, numCells, noiseLevel, templDSForProfile)
% DSAVEGenerateSNODataset
%   Generates a Sampling Noise Only (SNO) dataset.
%   This function uses a template dataset to compute the probability p that
%   a molecule from a certain gene should be picked each time a new UMI is
%   found. This is calculated from the counts value, i.e. counts/sum of all
%   counts. This is based on a multinomial distribution, and it will select
%   exactly the same number of UMIs that are in the original set using the
%   probabilities for the genes from the mean of the dataset sent in.
%   Important that this is really UMI counts in this case, not TPM!
%   Random multiplicative noise will be added (0 == no noise)
%   noiseLevel should be 0 or greater. A standard random normal distributed
%   noise multiplied by noiseLevel will be multiplied to the probabilities.
%   templDSForProfile - 
% Input:
%   templDs         The input dataset (cell population)
%   progrBarCtxt    (optional) Progress bar context.
%   numCells        (optional) Can be used to specify the number of cells. 
%                   Defaults to the number of cells in the input dataset.
%   noiseLevel      (optional) The noise level to add; defaults to 0 (no noise)
%   templDSForProfile (optional) can be used if you want to generate data 
%                   from a different cell type - defaults to templDs - the 
%                   genes in the datasets need to be synchronized
%
% Usage: ds = DSAVEGenerateSNODataset(ds)
%
% Johan Gustafsson, 2019-05-20
%
if nargin < 2
    progrBarCtxt = [];
end
if nargin < 3
    numCells = size(templDs.data,2);
end
if nargin < 4
    noiseLevel = 0;
end
if nargin < 5
   templDSForProfile = templDs; 
end
progbar = ProgrBar(['Generating SNO dataset from template ''' templDs.name ''''], progrBarCtxt);

%create empty dataset
ds = SCDataset;
ds.name = ['Sampling dataset from templ ' templDs.name ' for ' num2str(numCells) ' cells'];
ds.genes = templDs.genes;
numGenes = size(ds.genes,1);
ds.data = zeros(numGenes,numCells);%sparsify in the end
ds = ds.fillEmpties();


%generate probabilities
tmpTempl = TPM(templDSForProfile);
meanRefExpr = mean(tmpTempl.data,2);
%meanRefExpr(isnan(meanRefExpr)) = 0;
prob = meanRefExpr./sum(meanRefExpr,1);
%unsparse...
prob = full(prob);
%generate a vector with the sum of probabilities up to this gene (including this gene)
probSum = prob;%just create a vector of the right size
%probSum(1,1) = prob(1,1);%not really needed, but here for clarity
%for i = 2:numGenes
%    probSum(i,1) = probSum(i-1,1) + prob(i,1);
%end
probSum(1,1) = 0;%not really needed, but here for clarity
for i = 2:numGenes
    probSum(i,1) = probSum(i-1,1) + prob(i-1,1);
end

edges = [probSum;max(probSum(end,1),1)];%add right edge; make sure it is not smaller than the previous due to roundoff problems. 
%probSum(numGenes,1)
%prob(numGenes,1)


%create vector of number of UMIs
if (size(templDs.data,2) == numCells)
    UMIs = sum(templDs.data,1);
else
    UMIs = zeros(1, numCells);

    sumUMI = sum(templDs.data,1);
    ind = 1;
    cellsLeft = numCells;
    cellsInTemplate = size(templDs.data,2);
    while cellsLeft > 0
        cellsThisRound = min(cellsLeft, cellsInTemplate);
        sel = randsample(cellsInTemplate, cellsThisRound);
        UMIs(:,ind:ind+cellsThisRound-1) = sumUMI(1,sel);
        ind = ind + cellsThisRound;
        cellsLeft = cellsLeft - cellsThisRound;
    end

    UMIs = round(UMIs);%in case the dataset has been scaled or something
end

%generate cells
if noiseLevel == 0 % no noise, will go faster
    for i = 1:numCells
        if UMIs(1,i) > 0 %handle the fact that cells sometimes can contain 0 UMIs, which doesn't work with the code below. The data is automatically 0 in that case.
            progbar.Progress(i/numCells);
            %generate n random numbers between 0-1 and sort them
            r = rand([UMIs(1,i) 1]);
            ds.data(:,i) = histcounts(r,edges);
        end
    end
else %with noise
    numGenes = size(prob,1);
    for i = 1:numCells
        progbar.Progress(i/numCells);
        
        X = randn(numGenes, 1) * noiseLevel;
        noise = 2.^X;
        probTmp = prob .* noise;
        %scale probabilities back to the sum of 1
        probTmp = probTmp ./ sum(probTmp,1);
        
        %have to recalculate edges
        probSum(1,1) = 0;
        for j = 2:numGenes
            probSum(j,1) = probSum(j-1,1) + probTmp(j-1,1);
        end

        edges = [probSum;max(probSum(end,1),1)];%add right edge; need one more edge than num genes! 
        
        r = rand([UMIs(1,i) 1]);
        ds.data(:,i) = histcounts(r,edges);
    end
end

ds.data = sparse(ds.data);

progbar.Done();

end

