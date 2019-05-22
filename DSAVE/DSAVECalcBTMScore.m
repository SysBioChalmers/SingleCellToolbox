function [results,differenceCVs] = DSAVECalcBTMScore(origDs, templInfo, progrBarCtxt, skipAlignment, iterations, useLogTransform, logTPMAddon)
% DSAVECalcBTMScore
%   Aligns the dataset to the template, by removing genes, removing cells and
%   down-sampling the data to match the template as good as possible.
%   All datasets successfully aligned to the same template will have almost
%   identical sampling noise.
% Input:
%   origDs          The input dataset (cell population)
%   templInfo       The template to align to
%   progrBarCtxt    (optional) Progress bar context.
%   skipAlignment   (optional) If true, alignment is skipped
%   iterations      (optional) Number of iterations; defaults to 15
%   useLogTransform (optional) If true, calculation will be done on 
%                   log transformed data. Defaults to false.
%   logTPMAddon     (optional) When using log transformed data, this number
%                   is added to the data before transforming to avoid
%                   log(0).
% Output:
%   results         struct containing the score and other variation info
%   differenceCVs   difference CVs of the individual runs - can be used for 
%                   evaluating the number of iterations needed
% Usage: res = DSAVECalcBTMScore(ds, templInfo)
%
% Johan Gustafsson, 2019-05-20
%

if nargin < 3
    progrBarCtxt = ProgrBarContext;
end

if nargin < 4
    skipAlignment = false;
end

if nargin < 5
    iterations = 15;
end

if nargin < 6
    useLogTransform = false;
end

if nargin < 7
    logTPMAddon = 1;
end

progbar = ProgrBar(['Calculating DSAVE score for dataset  ''' origDs.name ''''], progrBarCtxt);


lbnonlog = 10;
ubnonlog = 1000;
numPoints = 1000;%use the same number of points to simplify for data allocation
meansLog = linspace(log10(lbnonlog), log10(ubnonlog), numPoints);
xes = 10 .^ meansLog;


numXVals = size(xes,2);

alignedCVs = zeros(iterations,numXVals);
samplingCVs = zeros(iterations,numXVals);
differenceCVs = zeros(iterations,numXVals);


for it = 1:iterations
    if (skipAlignment)
        aligned = origDs;
    else
        aligned = DSAVEAlignDataset(origDs, templInfo, progbar.GetSubContext(0.59/iterations));%use a different subset of cells in each loop iteration
    end
    
    SNO = DSAVEGenerateSNODataset(aligned, progbar.GetSubContext(0.40/iterations));
    if useLogTransform
        alignedGeneCVs = GetGeneCVsLog(aligned, templInfo, logTPMAddon);
        SNOGeneCVs = GetGeneCVsLog(SNO, templInfo, logTPMAddon);
    else
        alignedGeneCVs = GetGeneCVs(aligned, templInfo);
        SNOGeneCVs = GetGeneCVs(SNO, templInfo);
    end
    
    %throw away the most and least variable genes
    numGenes = size(aligned.data,1);
    numToDiscardUpper = round(templInfo.fractionUpperOutliers * numGenes);
    numToDiscardLower = round(templInfo.fractionLowerOutliers * numGenes);
    difference = alignedGeneCVs - SNOGeneCVs;
    [~,gi] = sort(difference);
    discard = [gi(1:numToDiscardLower); gi((end-numToDiscardUpper+1):end)];
    anyToDiscard = sum(discard,1) > 0;
    
    alignedGeneCVsRem = alignedGeneCVs;
    SNOGeneCVsRem = SNOGeneCVs;
    if anyToDiscard
        alignedGeneCVsRem(discard,:) = [];
        SNOGeneCVsRem(discard,:) = [];
    end
    
    alData = TPM(aligned.data);
    alData(discard,:) = [];
    SNOData = TPM(SNO.data);
    SNOData(discard,:) = [];
    
    
    [alCVs, alXes] = AverageIntoBins(alData, alignedGeneCVsRem, templInfo);
    [saCVs, saXes] = AverageIntoBins(SNOData, SNOGeneCVsRem, templInfo);

    %remove any NaNs (let linear interpolation fill in)
    
    %now, recalculate the x:es and cvs to the Xes in the template using linear interpolation;
    %otherwise it will be difficult to take the mean later
    
    %first remove any points with duplicate x; these will
    %otherwise mess up the linear interpolation.
    [alXes,ia,~] = unique(alXes);
    alCVs = alCVs(1,ia);
    [saXes,ia,~] = unique(saXes);
    saCVs = saCVs(1,ia);
    %then use linear interpolation
    %if the code fails on either of these two lines, that is because there are no genes that fits some bins
    %in that case, try using a template that discards fewer outliers
    alignedCVs(it,:) = interp1(alXes,alCVs,xes);
    samplingCVs(it,:) = interp1(saXes,saCVs,xes);
    
    %differenceCVs(it,:) = alignedCVs(it,:) - samplingCVs(it,:);
    differenceCVs(it,:) = alignedCVs(it,:) - samplingCVs(it,:);
    
end

%now, take the mean of the iterations from all three curves. 

results = struct();
results.alignedCVs = mean(alignedCVs,1);
results.samplingCVs = mean(samplingCVs,1);
results.differenceCVs = mean(differenceCVs,1);
results.DSAVEScore = mean(results.differenceCVs);%mean over all points ranging over different TPM
results.tpms = xes;

progbar.Done();

end


function logcv = GetGeneCVs(ds, templInfo)
    ds_red = TPM(ds);

    totset = ds_red.data;
    numGenes = size(totset, 1);

    avgRefExpr = mean(totset,2);
    %calc variances
    variances = var(totset, 0, 2);
    sd = sqrt(variances);
    cv_ = sd ./ (avgRefExpr + 0.05);%Coefficient of Variation = std/mean. Adding 0.05, a neglectably small number, to handle too lowly expressed genes
    logcv = log(cv_ + 1);%the + 1 says that no variance -> 0 value
end


function logcv = GetGeneCVsLog(ds, templInfo, logTPMAddon)
    ds_red = TPM(ds);
    ds_red.data = LogTrans(ds.data, 1, logTPMAddon);
    totset = ds_red.data;
    avgRefExpr = mean(totset,2);
    %calc variances
    variances = var(totset, 0, 2);
    sd = sqrt(variances);
    logcv = sd ./ (avgRefExpr + 0.01);%Coefficient of Variation = std/mean. Adding 0.05, a neglectably small number, to handle too lowly expressed genes
end



function [cv, meanGeneExpr] = AverageIntoBins(TPMData, logcv, templInfo)

avgRefExpr = mean(TPMData,2);


numBins = size(templInfo.binningInfo.means,2);

cv = zeros(1,numBins);
meanGeneExpr = templInfo.binningInfo.means;


for i = 1:numBins
    %select the genes within the expression range
    sel = avgRefExpr >= templInfo.binningInfo.lbs(1,i) & avgRefExpr <= templInfo.binningInfo.ubs(1,i);
    cv(1, i) = mean(logcv(sel));%y value in the graph
    meanGeneExpr(1,i) = 2^mean(log2(avgRefExpr(sel)+0.05)) - 0.05;%geometric-ish mean
end


end


