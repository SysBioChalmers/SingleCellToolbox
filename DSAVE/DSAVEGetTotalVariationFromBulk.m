function res = DSAVEGetTotalVariationFromBulk(s, pool4samples, upperBoundTPM, lowerBoundTPM)
% DSAVEGetTotalVariationFromBulk
%   Calculates the average pairwise total variation between a list of bulk samples, with a TPM
%   filtration of the genes.
% Input:
%   s               The samples to be investigated
%   pool4samples    If true, the algorithm compares the means of two groups
%                   of samples with 4 samples in each group. s is
%                   expected to have 8 samples in total for this option.
%   upperBoundTPM   All genes above this threshold are discarded.
%   lowerBoundTPM   All genes below this threshold are discarded.
%
% Usage: templInfo = DSAVEGetTotalVariationFromBulk(s, false, 10000000, 0);
%
% Johan Gustafsson, 2019-05-21
%

    means = mean(s.data,2);
    badGenes = means < lowerBoundTPM | means > upperBoundTPM;
    s = s.geneSubset(s.genes(~badGenes)); 
    numSamp = size(s.sampleIds,2);
    numGenes = size(s.genes,1);
    diffs = zeros(numGenes,numSamp*(numSamp-1)/2);
    index = 1;

    if ~pool4samples
        %loop though each pair of the samples
        for i = 1:numSamp-1
           for j = i+1:numSamp
               diffs(:,index) = log((s.data(:,i)+0.05)./(s.data(:,j)+0.05));
               index = index + 1;
           end
        end
    else
        if size(s.data,2) ~= 8
            error('4 on 4 only works with 8 samples');
        end
        %here the average of 4 samples is compared to the
        %average of another 4.
        %get all permutations of four
        ind = 1:numSamp;
        k = round(numSamp/2);
        combs = nchoosek(ind,k);
        numCombs = size(combs,1);
        diffs = zeros(numGenes,numCombs);
        for i = 1:numCombs
            a = s.data(:,combs(i,:));
            notind = 1:numSamp;
            notind(combs(i,:)) = [];
            b = s.data(:,notind);
            ma = mean(a,2);
            mb = mean(b,2);
            diffs(:,i) = log((ma+0.05)./(mb+0.05));
        end
    end
    res = mean(mean(abs(diffs),2),1);
end
