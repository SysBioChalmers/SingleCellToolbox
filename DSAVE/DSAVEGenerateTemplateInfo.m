function templInfo = DSAVEGenerateTemplateInfo(ds, datasetsForGenes, numCells, numUMIs, fractionUpperOutliers, fractionLowerOutliers, progrBarCtxt)
% DSAVEGenerateTemplateInfo
%   Generates a template that datasets can be aligned to.
% Input:
%   ds              The dataset (cell population) from which the UMI counts 
%                   per cell are taken
%   datasetsForGenes The genes in the template will be the intersection of
%                   the genes in the datasets of this horizontal cell array.
%   numCells        The number of cells to include in the template. 
%                   Defaults to the number of cells in the input dataset.
%   numUMIs         A random number of UMIs will be discarded from ds until
%                   the average number of UMIs per cell equals numUMIs
%   fractionUpperOutliers Sets the fraction of the most variable genes that
%                   should be discarded.
%   fractionLowerOutliers Sets the fraction of the least variable genes that
%                   should be discarded.
%   progrBarCtxt    (optional) Progress bar context.
%
% Usage: templ = DSAVEGenerateTemplateInfo(ds, {ds1, ds2, ds3}, 1000, 1000, 0.025, 0.025);
%
% Johan Gustafsson, 2019-05-20
%

if nargin < 7
    progrBarCtxt = [];
end
%% Binning info

%place the means evenly distributed in log scale
%we're interested in 10-1000; make sure we cover that range since the x:es
%are not always that accurate
lbnonlog = 5;
ubnonlog = 1300;
numPoints = 1000;
meansLog = linspace(log10(lbnonlog), log10(ubnonlog), numPoints);
templInfo.binningInfo.means = 10.^meansLog;

poolSize = 500;

%get mean values of the genes
gm = TPM(mean(ds.data,2));
gmLog = log10(gm);

templInfo.binningInfo.lbs = zeros(1, numPoints);
templInfo.binningInfo.ubs = zeros(1, numPoints);

progbar = ProgrBar('Generating template', progrBarCtxt);

for i = 1:numPoints
    progbar.Progress(i/numPoints * 0.9);
    %now increase the bin size until at least 100 genes are included
    step = 0.0005;
    lb = meansLog(1,i) - step;
    ub = meansLog(1,i) + step;
    for s = 1:1000
        sel = gmLog >= lb & gmLog <= ub;
        numGenes = sum(sel);
        if numGenes >= poolSize
            break;
        end
        %increase the bounds in a way that conserves the mean
        meanOfCaughtGenes = mean(gmLog(sel));
        if (meanOfCaughtGenes > meansLog(i))
            lb = lb - step;
        else
            ub = ub + step;
        end
    end
    templInfo.binningInfo.lbs(1,i) = 10^lb;
    templInfo.binningInfo.ubs(1,i) = 10^ub;
end


%% Geneset
dss = SynchronizeGenes(datasetsForGenes, [], true);
templInfo.geneSet = dss{1}.genes;

progbar.Progress(0.95);


%% UMIDistr
%first remove the genes that should not be included
subDs = ds.geneSubset(templInfo.geneSet);

%then select a subset of the cells without replacement
subDs = subDs.randSample(numCells);
sumUMIs = sum(subDs.data,1);

%then reduce the number of UMIs until we have the desired average:
%figure out how many umis to remove
totUMIs = sum(sumUMIs);
totTargetUMIs = numCells * numUMIs;
toRem = totUMIs - totTargetUMIs;
templInfo.UMIDistr = sumUMIs;%resulting UMIs per cell
if toRem < 0
    error('warning: Could not reach target UMIs per cell due to that the UMIs in the template dataset were to few');
elseif toRem > 0 %no point doing anything if the UMI count is already right
    %randomly select toRem number of UMIs to remove
    indToRem = randsample(totUMIs, toRem);
    edges = [0 (cumsum(templInfo.UMIDistr,2)+0.01)];%+0.01 changes the check on the upper bound from < to <= and lower bound from >= to >
    subtr = histcounts(indToRem,edges);
    templInfo.UMIDistr = templInfo.UMIDistr - subtr;%There is a risk to get cells with 0 UMIs, but we believe this risk is very small and can be neglected
end
templInfo.UMIDistr = sort(templInfo.UMIDistr);

%% Some other params
templInfo.fractionUpperOutliers = fractionUpperOutliers;
templInfo.fractionLowerOutliers = fractionLowerOutliers;

progbar.Done();

end
