function [genes, logCVDifference, pVals, SNOVariances, SNOCountsPerGene] = DSAVEGetGeneVariation(ds, lb, iterations, maxNumCells, progrBarCtxt)
% DSAVEGetGeneVariation
%   Calculates the DSAVE gene-wise BTM variation metric.
% Input:
%   ds              The dataset (cell population) to be investigated
%   lb              (optional) TPM lower bound, genes below this will not be
%                   investigated. Defaults to 10 TPM/CPM.
%   iterations      (optional) The number of SNO datasets to generate. 
%                   Recommended value is 100 if no p values are needed, 
%                   10,000 - 100,000 if p-values are of interest. Defaults
%                   to 100,000.
%   maxNumCells     (optional) ds is reduced to this number of cells if it
%                   contains more, to save computation time. Defaults to 10,000. 
%   progrBarCtxt    (optional) Progress bar context.
%
% Output:
%   genes           The genes that passed the lower bound filter. The order
%                   of the genes here is the same as the order in
%                   logCVDifference and pVals.
%   logCVDifference The BTM variation for each gene.
%   pVals           The p values for each gene.
%   SNOVariances    The variances for all iterations. As an optimization,
%                   only a single row is stored for each total count value
%                   per gene in the dataset. Thus, this does not match the
%                   genes in the genes output variable. The size of this
%                   variable can be very large. A standard user normally do
%                   not need this output variable.
%   SNOCountsPerGene The total counts per gene used for each row in 
%                   SNOVariances. Each gene need to be mapped to these
%                   counts before being able to use SNOVariances. A standard 
%                   user normally do not need this output variable.
%
% Usage: templ = DSAVEGetGeneVariation(ds, 1, 10000, 1000);
%
% Johan Gustafsson, 2019-05-20
%
if nargin < 5
    progrBarCtxt = [];
end
if nargin < 4
    maxNumCells = 2000;%cut down to 2000 cells to get somewhat reasonable computation times
end
if nargin < 3
    iterations = 100000;
end
if nargin < 2
    lb = 10;
end

progbar = ProgrBar(['Calculating gene-wise variation for dataset ''' ds.name ''''], progrBarCtxt);


numCells = size(ds.data,2);
if (numCells > maxNumCells)
    ds = ds.randSample(maxNumCells);
    numCells = maxNumCells;%keep this line if we start using numCells below later
end

dstpm = TPM(mean(ds.data, 2));

sel = dstpm >= lb & dstpm ~= 0 & sum(ds.data,2) ~= 1; %skip all below threshold, and all with 0 or 1 counts, doesn't make sense to look at those 

ds = ds.geneSubset(sel);%filter genes to avoid problems with division by zero, etc, and to save compilation time. No point in looking at too lowly expressed ones anyway, they will not become significant.
numGenes = size(ds.data,1);

%generate SNO TPMs
%SNOTPMs = (10.^(-.3:0.005:4)).';
%convert the TPMs into counts
SNOUMIsPerCell = sum(ds.data,1);
%sumUMIsPerCell = sum(SNOUMIsPerCell,2);

SNOCountsPerGene = full(unique(sum(ds.data,2)));
%skip 0 and 1
SNOCountsPerGene(SNOCountsPerGene < 2) = [];

%countsPerGenePre = round(SNOTPMs.*sumUMIsPerCell./10^6);
%SNOCountsPerGene = unique(countsPerGenePre); %some of the low TPMs are sometimes the same count

numCountVals = size(SNOCountsPerGene,1);


SNOLogCVS = zeros(numCountVals, iterations);
SNOVariances = zeros(numCountVals, iterations);
[logCVDS,varianceDS] = GetVarAndLogCV(ds.data, sum(ds.data,1));

%precalc things for generating SNO datasets
%generate probabilities
prob = SNOUMIsPerCell./sum(SNOUMIsPerCell,2);
prob = full(prob);

%generate a vector with the sum of probabilities up to this cell (including this cell)
probSum = prob;%just create a vector of the right size
probSum(1,1) = 0;%not really needed, but here for clarity
for ii = 2:numCells
    probSum(1, ii) = probSum(1,ii-1) + prob(1, ii-1);%can probably speed this up with cumsum
end

SNOEdges = [probSum max(probSum(1,end),1)];%add right edge; make sure it is not smaller than the previous due to roundoff problems. 
SNOEdges = SNOEdges.';
SNOdata = zeros(numCountVals,numCells);%preallocate data matrix that will be reused and overwritten in each iteration, to avoid copying of data

for it = 1:iterations
    SNOdata = GenSampDs(SNOdata, SNOCountsPerGene, SNOEdges);
    [SNOLogCVS(:,it),SNOVariances(:,it)] = GetVarAndLogCV(SNOdata, SNOUMIsPerCell);%so, use the UMIs per cell from the original dataset when TPM:ing
    progbar.Progress(it/iterations);
end

genes = ds.genes;

%Calculate p value with a non-parametric method:

countsPerGene = full(sum(ds.data,2));
SNOCountsPerGene = full(SNOCountsPerGene);

%map counts to SNO counts
indices = arrayfun( @(x)( find(SNOCountsPerGene==x) ), countsPerGene);%should only find one index per row

sz = size(SNOVariances, 2);

pVals = sum(SNOVariances(indices,:) >= varianceDS,2) ./ sz;


%so subtraction of CVs
logCVSNOm = mean(SNOLogCVS,2);

logCVDifference = logCVDS - logCVSNOm(indices,:);



progbar.Done();

end

function [logCV,variances] = GetVarAndLogCV(dsdata, UMIsPerCell)
    dsdata = dsdata .* 10^6 ./ UMIsPerCell;
    avgRefExpr = mean(dsdata,2);
    variances = var(dsdata, 0, 2);
    sd = sqrt(variances);
    cv_ = sd ./ (avgRefExpr + 0.05);%Coefficient of Variation = std/mean. Adding 0.05, a neglectably small number, to handle too lowly expressed genes
    logCV = log(cv_ + 1);%the + 1 says that no variance -> 0 value
end
 
%Generates a SNO dataset where gene expression is preserved instead of
%counts per cell. Note that the data structure is reused, need to be
%allocated outside this function and used according to data =
%data = GenSampDs(data, ...)
function data = GenSampDs(data, countsPerGene, edges)
    numGenes = size(data,1);
    %generate genes
    for i = 1:numGenes
        if countsPerGene(i,1) > 0 %handle the fact that cells sometimes can contain 0 UMIs, which doesn't work with the code below. The data is automatically 0 in that case.
            %generate n random numbers between 0-1 and sort them
            r = rand([countsPerGene(i,1) 1]);
            a = histcounts(r,edges);
            data(i,:) = a.';
        end
    end
end


