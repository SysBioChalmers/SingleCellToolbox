function s = PoolSC2Samples(scs)
% Takes a cell array of SCDatasets, pools the cells of each of them and
% puts the results in a Samples object. The genes are expected to be
% synchronized.
%
% Input:
%   scs  Cell array of SCDatasets
% Usage: samp = PoolSC2Samples({ds1,ds2,ds3, ...});
%
% Johan Gustafsson, 2019-09-15
%
    numSamp = size(scs,2);
    numGenes = size(scs{1,1}.genes, 1);
    s = Samples;
    s.data = zeros(numGenes, numSamp);
    s.genes = scs{1,1}.genes;
    s.sampleIds = cell(1,numSamp);
    
    for i = 1:numSamp
        s.sampleIds{1,i} = scs{1,i}.name;
        s.data(:,i) = mean(TPM(scs{1,i}.data),2);
    end
end


