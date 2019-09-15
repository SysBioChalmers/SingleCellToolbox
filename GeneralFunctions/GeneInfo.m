classdef GeneInfo
   properties
      %these only contain genes from the normal chromosomes 1-22, X and Y
      ensGeneId %ensembl stable id without version, vertical cell array
      geneName%gene name, i.e. CD276 etc, vertical cell array
      chromosome%chromosome, vertical cell array
      geneStart%gene start
      geneEnd%gene end
      meanTranscriptLength%average transcript length for all splice variants
      
      %copies of the above where all other chromosomes exist as well, so
      %many gene names are duplicated even though their ensembl stable ids are
      %unique
      %this gives a better coverage when converting from ensembl stable id
      %to gene id, but doesn't help in the other direction
      ensGeneIdFull %ensembl stable id without version, vertical cell array
      geneNameFull%gene name, i.e. CD276 etc, vertical cell array
      chromosomeFull%chromosome, vertical cell array
      geneStartFull%gene start
      geneEndFull%gene end
      meanTranscriptLengthFull%average transcript length for all splice variants
   end
   methods
      function gn = Ens2GeneName(this, A, clearGenesWithoutMatch)%A is expected to be a vertical cell array 
          if (nargin < 3)
              clearGenesWithoutMatch = false;
          end
          A = RemoveVersionOnEns(this, A);
          [~, ia, ib] = intersect(A, this.ensGeneIdFull);
          [r,~] = size(A);
          gn = cell(r,1);
          gn(ia) = this.geneNameFull(ib);
          if ~clearGenesWithoutMatch
            empties = cellfun(@isempty,gn);
            gn(empties) = A(empties); 
          end
      end 
      function gn = GeneName2Ens(this, A)%A is expected to be a vertical cell array 
          %not done
          [~, ia, ib] = intersect(A, this.geneName);
          [r,~] = size(A);
          gn = cell(r,1);
          gn(ia) = this.ensGeneId(ib);
          empties = cellfun(@isempty,gn);
          gn(empties) = A(empties); 
      end 
      %data for genes not recognized will be left untouched
      function res = NormalizeByTranscriptLength(this, ds)
          res = ds;
          [~, ia, ib] = intersect(ds.genes, this.geneName);
          avg = mean(this.meanTranscriptLength(ib));
          res.data(ia,:) = (res.data(ia,:) ./ this.meanTranscriptLength(ib)) .* avg;
      end 
      function gn = RemoveVersionOnEns(this, A)
        tokens = regexp(A, '([^.])*\.?[0-9]*', 'tokens');
        ensGenesTmp = [tokens{:}];
        gn = [ensGenesTmp{:}].';
      end
   end
end