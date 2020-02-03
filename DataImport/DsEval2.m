classdef DsEval2
    % DsEval2
    %   Reads Eval2 data from file into an SCDataset class, GSE75790, mouse
    %   data. Contains ERCC genes. The genes are in ensembl format.
    %   The dataset contains data from many techniques, in two batches per
    %   technique. Note that the downloaded text file has been modified - I
    %   added a tab in the upper left corner to simplify the reading.
    %
    % Johan Gustafsson, 2020-01-24
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsEval2.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading Eval2...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/eval2.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsEval2.import('../../ImportableData/Eval2/GSE75790_ziegenhain_complete_data.txt');
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
        function  ds = import(path)
            % import
            %   Imports the data
            % Input:
            %   path       Path to the data file
            %
            % Usage: ds = DsEval2.import('../../ImportableData/Eval2/GSE75790_ziegenhain_complete_data.txt');
            %
            
            ds = SCDataset;
            ds.name = 'eval2';
            L = importdata(path,'\t');
            ds.data = L.data;
            [m,n] = size(ds.data);
            ds.cellIds = L.textdata(1, 2:(n+1));
            ds.genes = L.textdata(2:(m+1), 1);
            %extract method and put in 
            ds.sampleIds = arrayfun(@DsEval2.extr, ds.cellIds, 'UniformOutput', false);
            ds = ds.fillEmpties();
        end
        
        function res = extr(x)
            a = split(x,'_');
            res = a{1};
        end
        
    end
end

