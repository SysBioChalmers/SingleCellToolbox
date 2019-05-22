classdef Samples
    % Samples
    %   Class representing a several samples.
    %   For passing values by reference, follow the following pattern:
    %   s = s.fillEmpties(); % will avoid copying of data
    %
    % Johan Gustafsson, 2019-05-21
    %
    properties
        data %columns are samples, rows are genes
        genes % one column with all gene names. Use name convention "GAPDH" etc and not "ENS..." or "NM..."
        sampleIds % one row with sample ids
    end
    methods
        function s = sampleSubset(this, l) %logical vector
            % sampleSubset
            %   Extracts a subset of samples from a dataset into a new Samples
            %   object
            % Input:
            %   l               Logical vector or a vector of indices
            %
            % Usage:
            %   s2 = s.sampleSubset([1,3,8:23]);
            %
            s = Samples;
            s.data = this.data(:, l);
            s.genes = this.genes;
            s.sampleIds = this.sampleIds(1, l);
        end
        
        function [s,ia] = geneSubset(this, genesToKeep, sortGenes)
            % geneSubset
            %   Creates a new dataset only containing certain genes.
            % Input:
            %   genesToKeep     The genes to keep. Should be a vertical
            %                   cell array, logical array or array of indices
            %   sortGenes       If true, sorts genes, otherwise keeps the
            %                   order of the genes.
            %
            % Usage:
            %   s2 = s.geneSubset({'GAPDH'; 'RPS10'});
            %
            if nargin < 3
                sortGenes = 0;
            end
            s = this;
            if isnumeric(genesToKeep) || islogical(genesToKeep)
                s.genes = this.genes(genesToKeep);
                s.data = s.data(genesToKeep,:);
                if sortGenes
                    error('Sort genes are not supported for input of logical array or indices');
                end
            else
                if sortGenes
                    [s.genes,ia,~] = intersect(this.genes, genesToKeep, 'sorted');
                else
                    [s.genes,ia,~] = intersect(this.genes, genesToKeep, 'stable');
                end
                s.data = s.data(ia,:);
            end
            
        end
        
        function sRes = innerJoin(this, s)
            % innerJoin
            %   Joins two samples objects, discarding any genes that do not 
            %   exist in both.
            % Input:
            %   s              The dataset to join with
            %
            % Usage: s3 = s.innerJoin(s2);
            %
            sRes = Samples;
            [sRes.genes, ia, ib] = intersect(this.genes, s.genes);
            sRes.data = [this.data(ia,:) s.data(ib,:)];
            sRes.sampleIds = [this.sampleIds s.sampleIds];
        end
        
        function sRes = fullOuterJoin(this, s)
            % fullOuterJoin
            %   Joins two samples objects, keeping any genes that do not 
            %   exist in both objects. Missing data is imputed with 0.
            % Input:
            %   s              The dataset to join with
            %
            % Usage: s3 = s.fullOuterJoin(s2);
            %
            sRes = Samples;
            sRes.genes = union(this.genes, s.genes);
            %join the data matrices and adapt to the new gene set
            [~,ires,ids] = intersect(sRes.genes,s.genes);
            [numGenes,~] = size(sRes.genes);
            [~, numSamples] = size(s.sampleIds);
            sData = zeros(numGenes, numSamples);
            sData(ires,:) = s.data(ids,:);
            [~,ires,ithis] = intersect(sRes.genes,this.genes);
            [~, numSamples] = size(this.sampleIds);
            thisData = zeros(numGenes, numSamples);
            thisData(ires,:) = this.data(ithis,:);
            sRes.data = [thisData sData];
            sRes.sampleIds = [this.sampleIds s.sampleIds];
        end
        
        %We could also implement left join, which would keep the genes of
        %this and fill in with zeros for ds if missing, and discard the
        %extra genes from the ds set. If we need it...
        
        function s = fillEmpties(this)
            % fillEmpties
            %   Fills in all empty member variables with default values
            % Input:
            %   filename        Filename
            %
            % Usage: ds = ds.fillEmpties();
            %            
            s = this;
            [r,c] = size(s.data);
            if isempty(s.sampleIds)
                formatString = strcat('%d');
                s.sampleIds = sprintfc(formatString,(1:size(s.data,2)));
            end
        end
        
        function writeToTextFile(this, filename)
            % writeToTextFile
            %   Saves the data table, including sample ids and genes, to
            %   file.
            % Input:
            %   filename        Filename
            %
            % Usage: ds.writeToTextFile('test.tsv');
            %            
            fileID = fopen(filename,'w');
            fprintf(fileID,'\t%s', this.sampleIds{:,:});
            fprintf(fileID,'\n');
            rows = size(this.data, 1);
            for i = 1:rows
                fprintf(fileID,'%s', this.genes{i});
                fprintf(fileID,'\t%f', this.data(i,:));
                fprintf(fileID,'\n');
            end
            fclose(fileID);
        end
    end
end