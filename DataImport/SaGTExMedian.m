classdef SaGTExMedian
    % SaGTExMedian
    %   Reads the GTEx median tissue RNA-Seq bulk expression.
    %   It could be that the file to be read has been modified, and I removed
    %   the two top lines since they just caused problems, not sure.
    %
    % Johan Gustafsson, 2019-05-22
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the samples. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: s = SaGTExMedian.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading GTEx median...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/GTExMed.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = SaGTExMedian.import('../../ImportableData/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct');
                    save(filename, 'v');
                else
                    a = load(filename);
                    v = a.v;
                end
                DsHelper.restoreDir(prevDir);
            end
            ret = v;
        end
    end
    
    methods(Static, Access = private)
        function s = import(path)
            % import
            %   Imports the data
            % Input:
            %   directoryPath       Path to the 10x files. No slash at the end.
            %
            % Usage: s = SaGTExMedian.import('../../ImportableData/PBMCCD4TCellsMemory/filtered_matrices_mex/hg19');
            %
            
            %path = '../../ImportableData/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct';
            %path = 'C:/Work/MatlabCode/components/SCLib/ImportableData/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct';
            s = Samples;
            
            L = importdata(path,'\t');
            s.data = L.data;
            s.sampleIds = L.textdata(1, 3:end);
            s.genes = L.textdata(2:end, 2);
            s = TPM(s);
        end
    end
end
