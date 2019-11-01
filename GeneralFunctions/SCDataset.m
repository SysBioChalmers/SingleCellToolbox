classdef SCDataset
    % SCDataset
    %   Class representing a single-cell dataset.
    %   For passing values by reference, follow the following pattern:
    %   ds = ds.setSubCellType(...); % will avoid copying of data
    %
    % Johan Gustafsson, 2019-05-20
    %
    properties
        name % The name of the dataset
        data % Genes as rows, samples as columns
        genes % One column with all gene names. Use name convention "GAPDH" etc and not "ENS..." or "NM..."
        cellIds % One row of all cell ids
        sampleIds % One row of all tumors/patients/some other group.
        cellType % One row describing cell type. Shall match the values in the Celltype class.
        subCellType % One row which can be used for dividing the cells within a celltype into more categories
        extraCellInfo % One row (cell array of char arrays) that can be used to describe the cell in more detail
    end
    methods
        function ds = cellSubset(this, l)
            % cellSubset
            %   Extracts a cell subset of a dataset into a new dataset
            % Input:
            %   l               l is a logical vector or a vector of indices
            %
            % Usage: 
            %   ds2 = ds.cellSubset([1,3,8:23]); % gets specific cells
            %   ds2 = ds.cellSubset(ds.cellType == Celltype.BCell); % gets all B cells
            %
            ds = SCDataset;
            ds.name = strcat(this.name, '_subset');
            ds.data = this.data(:, l);
            ds.genes = this.genes;
            ds.cellIds = this.cellIds(1, l);
            ds.sampleIds = this.sampleIds(1, l);
            if (~isempty(this.cellType))
                ds.cellType = this.cellType(1, l);
            end
            if (~isempty(this.subCellType))
                ds.subCellType = this.subCellType(1, l);
            end
            if (~isempty(this.extraCellInfo))
                ds.extraCellInfo = this.extraCellInfo(1, l);
            end
        end
        
        function ds = randSample(this, numCells)
            % randSample
            %   Randomly selects numCells cells from the dataset, without
            %   replacement.
            % Input:
            %   numCells        The number of cells to extract
            %
            % Usage: 
            %   ds2 = ds.randSample(1000);
            %
            indToKeep = randsample(size(this.data,2),numCells);
            ds = this.cellSubset(indToKeep);
        end
        
        function [ds,ia] = geneSubset(this, genesToKeep) %genesToKeep should be a vertical cell array, logical array or array of indices
            % geneSubset
            %   Creates a new dataset only containing certain genes.
            % Input:
            %   genesToKeep     The genes to keep. Should be a vertical 
            %                   cell array, logical array or array of indices
            %
            % Usage: 
            %   ds2 = ds.geneSubset({'GAPDH'; 'RPS10'});
            %
            ds = this;
            if isnumeric(genesToKeep) || islogical(genesToKeep)
                ds.genes = this.genes(genesToKeep);
                ds.data = ds.data(genesToKeep,:);
            else
                [ds.genes,ia,~] = intersect(this.genes, genesToKeep, 'stable');
                ds.data = ds.data(ia,:);
            end
        end
        
        function ds = setCellType(this, cellIds, values)
            % setCellType
            %   Shortcut for setting cell type
            % Input:
            %   cellIds         Cell ids of cells to update 
            %   values          the cell types to set
            %
            % Usage: 
            %   ds = ds.setCellType(cellIds, values);
            %
            ds = this;
            [~,ia,ib] = intersect(ds.cellIds, cellIds);
            ds.cellType(ia) = values(ib);
        end
        
        function ds = setSubCellType(this, cellIds, values)
            % setSubCellType
            %   Shortcut for setting sub cell type
            % Input:
            %   cellIds         Cell ids of cells to update 
            %   values          the sub cell types to set
            %
            % Usage: 
            %   ds = ds.setSubCellType(cellIds, values);
            %
            ds = this;
            [~,ia,ib] = intersect(ds.cellIds, cellIds);
            ds.subCellType(ia) = values(ib);
        end
        
        
        function dsRes = innerJoin(this, ds)
            % innerJoin
            %   Joins two datasets, discarding any genes that do not exist in
            %   both datasets.
            % Input:
            %   ds              The dataset to join with
            %
            % Usage: ds3 = ds.innerJoin(ds2);
            %
            dsRes = SCDataset;
            dsRes.name = strcat('inner join (', this.name,', ',ds.name, ')');
            [dsRes.genes, ia, ib] = intersect(this.genes, ds.genes);
            dsRes.data = [this.data(ia,:) ds.data(ib,:)];
            dsRes.cellIds = [this.cellIds ds.cellIds];
            dsRes.sampleIds = [this.sampleIds ds.sampleIds];
            dsRes.cellType = [this.cellType ds.cellType];
            dsRes.subCellType = [this.subCellType ds.subCellType];
            dsRes.extraCellInfo = [this.extraCellInfo ds.extraCellInfo];
        end
        
        function dsRes = fullOuterJoin(this, ds)
            % fullOuterJoin
            %   Joins two datasets, keeping any genes that do not exist in
            %   both datasets. Missing data is imputed with 0.
            % Input:
            %   ds              The dataset to join with
            %
            % Usage: ds3 = ds.fullOuterJoin(ds2);
            %
            dsRes = SCDataset;
            dsRes.name = strcat('full outer join (', this.name,', ',ds.name, ')');
            dsRes.genes = union(this.genes, ds.genes);
            %join the data matrices and adapt to the new gene set
            [~,ires,ids] = intersect(dsRes.genes,ds.genes);
            [numGenes,~] = size(dsRes.genes);
            [~, numCells] = size(ds.cellIds);
            dsData = zeros(numGenes, numCells);
            dsData(ires,:) = ds.data(ids,:);
            [~,ires,ithis] = intersect(dsRes.genes,this.genes);
            [~, numCells] = size(this.cellIds);
            thisData = zeros(numGenes, numCells);
            thisData(ires,:) = this.data(ithis,:);
            dsRes.data = [thisData dsData];
            
            dsRes.cellIds = [this.cellIds ds.cellIds];
            dsRes.sampleIds = [this.sampleIds ds.sampleIds];
            dsRes.cellType = [this.cellType ds.cellType];
            dsRes.subCellType = [this.subCellType ds.subCellType];
            dsRes.extraCellInfo = [this.extraCellInfo ds.extraCellInfo];
        end
        
        %We could also implement left join, which would keep the genes of
        %this and fill in with zeros for ds if missing, and discard the
        %extra genes from the ds set. If we need it...
        
        function saveDataTable(this, filename)
            % saveDataTable
            %   Saves the data table, including cell ids and genes, to
            %   file.
            % Input:
            %   filename        Filename
            %
            % Usage: ds.saveDataTable('test.tsv');
            %            
            t = array2table(this.data);
            t.Properties.RowNames = this.genes;
            t.Properties.VariableNames = this.cellIds;
            writetable(t,filename,'Delimiter','\t','WriteRowNames',true);
        end
        
        function saveInfoTable(this, filename)
            % saveInfoTable
            %   Saves the cell types and sample ids to file
            % Input:
            %   filename        Filename
            %
            % Usage: ds.saveInfoTable('testinfo.tsv');
            %            
            t = table(this.cellType.', CelltypeId2CelltypeName(this.cellType).', this.sampleIds.','RowNames', this.cellIds.','VariableNames',{'celltype', 'celltype_text', 'sample_id'});
            writetable(t,filename,'Delimiter','\t','WriteRowNames',true);
        end
        
        function ds = fillEmpties(this)
            % fillEmpties
            %   Fills in all empty member variables with default values
            % Input:
            %   filename        Filename
            %
            % Usage: ds = ds.fillEmpties();
            %            
            ds = this;
            [r,c] = size(ds.data);
            if isempty(ds.cellIds)
                formatString = strcat(ds.name, '_%d');
                ds.cellIds = sprintfc(formatString,(1:size(ds.data,2)));
            end
            if isempty(ds.sampleIds)
                [ds.sampleIds{1:c}] = deal('Unknown');
            end
            if isempty(ds.cellType)
                ds.cellType = ones(1,c);%one happens to be Celltype.Unknown.
            end
            if isempty(ds.subCellType)
                ds.subCellType = zeros(1,c);
            end
            if isempty(ds.extraCellInfo)
                [ds.extraCellInfo{1:c}] = deal('');
            end
            
        end
        
        function samples = splitIntoRandomGroupSamples(this, groupSize)
            % splitIntoRandomGroupSamples
            %   Creates a samples object with random cell pools of size
            %   groupSize, without reusing cells. The algorithm will create 
            %   as many samples as possible from the dataset.
            % Input:
            %   groupSize        The size (number of cells) of the groups.
            %
            % Usage: samples = ds.splitIntoRandomGroupSamples(1000);
            %            
            totNumCells = size(this.cellIds,2);
            indices = randperm(totNumCells);
            
            numSamp = floor(totNumCells/groupSize);
            numGenes = size(this.genes,1);
            
            samples = Samples();
            samples.genes = this.genes;
            samples.sampleIds = cell(1,numSamp);
            samples.data = zeros(numGenes,numSamp);
            
            for i = 1:numSamp
                st = 1 + (i-1)*numSamp;
                en = i*numSamp;
                ind = indices(st:en);
                subds = this.cellSubset(ind);
                samples.sampleIds{1,i} = strcat(this.name,'_',num2str(i));
                samples.data(:,i) = mean(subds.data,2);
            end
        end
		
        function samples = poolToSingleSample(this)
            % getAsSingleSample
            %   Creates a samples object one sample, represented by the sum of all
            %   counts in the single cells. The pool is not converted to TPM by default.
            % Input:
            %
            % Usage: samples = ds.poolToSingleSample();
            %            

            samples = Samples();
            samples.genes = this.genes;
            samples.sampleIds = {"Pool"};
            samples.data = mean(this.data,2);
        end
    end
end