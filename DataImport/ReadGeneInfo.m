%Reads gene info from file
function gi = ReadGeneInfo(path)
gi = GeneInfo;
warning('OFF', 'MATLAB:table:ModifiedVarnames')
T = readtable(path, 'delimiter','\t');

ensGeneId = table2cell(T(:,1));
geneName = table2cell(T(:,7));
chromosome = table2cell(T(:,3));
geneStart = table2array(T(:,4));
geneEnd = table2array(T(:,5));
transcrLength = table2array(T(:,6));

%filter out the unique genes (they are duplicated for each transcript of the gene)
[gi.ensGeneIdFull, ia, ic] = unique(ensGeneId, 'stable');
gi.geneNameFull = geneName(ia);
gi.chromosomeFull = chromosome(ia);
gi.geneStartFull = geneStart(ia);
gi.geneEndFull = geneEnd(ia);

[a,~] = size(ia);
[c,~] = size(ic);
gi.meanTranscriptLengthFull = zeros(a,1);
%this takes forever, change it and use the fact that all identical ids come
%in after one another
lastId = '';
currIndex = 0;
sumTranscrLen = 0;
numSplices = 0;
for i = 1:c
    %check if we have reached a new gene or reached the end
    if ~strcmp(lastId, ensGeneId(i))
        %first store the collected average transcription lengths
        if currIndex ~= 0
            gi.meanTranscriptLengthFull(currIndex) = sumTranscrLen/numSplices;
        end
        %then restore all and step to the next 
        lastId = ensGeneId(i);
        currIndex = currIndex + 1;
        sumTranscrLen = 0;
        numSplices = 0;
    end
    sumTranscrLen = sumTranscrLen + transcrLength(i);
    numSplices = numSplices + 1;
end
%add the last element
gi.meanTranscriptLengthFull(currIndex) = sumTranscrLen/numSplices;

%now create the dataset part where all genes from the non-common chromosome
%set have been removed
commonChromosomes = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13'; ...
                     '14';'15';'16';'17';'18';'19';'20';'21';'22';'X';'Y'};
lia = ismember(gi.chromosomeFull, commonChromosomes);
gi.ensGeneId = gi.ensGeneIdFull(lia);
gi.geneName = gi.geneNameFull(lia);
gi.chromosome = gi.chromosomeFull(lia);
gi.geneStart = gi.geneStartFull(lia);
gi.geneEnd = gi.geneEndFull(lia);
gi.meanTranscriptLength = gi.meanTranscriptLengthFull(lia);
