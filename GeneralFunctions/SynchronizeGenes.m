function [res1,res2] = SynchronizeGenes(obj1, obj2, discardGenes)
% SynchronizeGenes
%   Makes sure the datasets/profiles/samples have the same genes in the same
%   order. Genes that don't exist in all are set to zero or discarded,
%   depending on the discardGenes flag
% Input:
%   obj1            First object or object list, horizontal cell array
%   ds2             Second objectTrue if transforming to log scale, false if from log
%                   scale.
%   numToAdd        (optional) The number to add to all data before log.
%                   Defaults to 1.
% Usage: Two possibilities:
%   1. [res1,res2] = SynchronizeGenes(ds1, ds2, discardGenes)
%   2. [res1,~] = SynchronizeGenes(dslist, [], discardGenes) - where dslist 
%   is a horizontal cell array of objects to be synchronized
%
% Johan Gustafsson, 2019-05-21
%

if (iscell(obj1))
    %obj 2 is expected to be empty
    if discardGenes
        %find intersection of all datasets; then synchronize on that
        isGenes = obj1{1,1}.genes;
        for i = 2:size(obj1,2)
            isGenes = intersect(isGenes, obj1{1,i}.genes);
        end

        %synchronize
        res1 = obj1;
        for i = 1:size(obj1,2)
            [res1{1,i}.genes, i1, ~] = intersect(obj1{1,i}.genes, isGenes);
            res1{1,i}.data = obj1{1,i}.data(i1,:);
        end
    else
        disp('not impl!!!!');
    end    
else
    res1 = obj1;
    res2 = obj2;
    if discardGenes
        [res1.genes, i1, i2] = intersect(obj1.genes, obj2.genes);
        res2.genes = res1.genes;
        res1.data = obj1.data(i1,:);
        res2.data = obj2.data(i2,:);
    else
        res1.genes = union(obj1.genes, obj2.genes);
        res2.genes = res1.genes;
        %join the data matrices and adapt to the new gene set
        numGenes = size(res1.genes,1);
        [~,ires,iobj] = intersect(res1.genes,obj1.genes);
        numCells = size(obj1.data,2);
        res1.data = zeros(numGenes, numCells);
        res1.data(ires,:) = obj1.data(iobj,:);

        [~,ires,iobj] = intersect(res2.genes,obj2.genes);
        numCells = size(obj2.data,2);
        res2.data = zeros(numGenes, numCells);
        res2.data(ires,:) = obj2.data(iobj,:);
    end
end

