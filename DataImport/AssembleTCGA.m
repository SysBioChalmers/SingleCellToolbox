%This file loads the melanoma TCGA files into a Samples structure
%note that the gene ids are in ensembl format, since conversion creates
%duplicates such as 'Y_RNA'
%Run:
%AssembleTCGA('C:/Work/MatlabCode/projects/AgileCancerTreatment/data/TCGAFiles/TCGA-LUSC', 'C:/Work/MatlabCode/components/SingleCellToolbox/TempData/tcga_lusc.txt', GeneInfo.get())

%AssembleTCGA('C:\Work\MatlabCode\projects\AgileCancerTreatment\data\TCGADownload\TCGA-LUSC\harmonized\Transcriptome_Profiling\Gene_Expression_Quantification')
function AssembleTCGA(path, outPath, geneInfo)

    myFiles = dir(fullfile(path,'*.txt')); %gets all wav files in struct
    remove = 0;
    for k = 1:length(myFiles)
        baseFileName = myFiles(k).name;
        fullFileName = fullfile(path, baseFileName);
        %fprintf(1, 'Now reading %s\n', fullFileName);
        str = sprintf('Reading file %d of %d', k, length(myFiles));
        strRem = '';
        for rem = 1:remove
            fprintf('\b');
            strRem = strcat(strRem,'\b');
        end
        fprintf('%s', str);
        remove = length(str);
        
        %for each file, read into an own samples structure and make a full
        %outer join with the empty one
        %this will be slow, but will be done only once... The problem is
        %that I don't know if the gene column will always be the same...
        tmpSamp = Samples();
        f = readtable(fullFileName, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
        tmpSamp.data = table2array(f(:,2));
        %remove '.FPKM.txt' from the file name to create a sample id
        tokens = regexp(baseFileName, '(.*)\.FPKM.txt', 'tokens');
        tmpSamp.sampleIds = [tokens{:}];
        
        ensGenesWithVer = table2cell(f(:, 1));
        %first we need to remove the version from the ensembl id. Then we
        %need to convert it to gene id.
        tokens = regexp(ensGenesWithVer, '([^.])*\.?[0-9]*', 'tokens');
        ensGenesTmp = [tokens{:}];
        ensGenes = [ensGenesTmp{:}].';
        tmpSamp.genes = geneInfo.Ens2GeneName(ensGenes);
        
        if k==1
            tcga_samp = tmpSamp;
        else
            tcga_samp = tcga_samp.fullOuterJoin(tmpSamp);
        end
    end
    
    tcga_samp.writeToTextFile(outPath);
end

