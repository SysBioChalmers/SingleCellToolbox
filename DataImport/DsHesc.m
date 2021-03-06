classdef DsHesc
    % DsHesc
    %   Reads hESC (human embryonic stem cell) data from file into an SCDataset class
    %   Not UMI based, few cells (in total 213 + 247 cells). This dataset
    %   is currently not tested at all! GSE64016
    %
    % Johan Gustafsson, 2020-01-22
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsHesc.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading human embryonic stem cells...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/hesc.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsHesc.import('../../ImportableData/hESC/GSE64016_H1andFUCCI_normalized_EC.csv');
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
            % Usage: ds = DsHesc.import('../../ImportableData/hESC/GSE64016_H1andFUCCI_normalized_EC.csv');
            %
            
            ds = SCDataset;
            ds.name = 'hesc';
            L = importdata(path,',');
            ds.data = L.data;
            [m,n] = size(ds.data);
            ds.cellIds = L.textdata(1, 2:(n+1));
            ds.genes = L.textdata(2:(m+1), 1);
           
            ds.cellType = repmat(Celltype.hESC,1,length(ds.cellIds));
            ds = ds.fillEmpties();
            ds = TPM(ds);
        end
    end        
end

