classdef DsMel
    % DsMel
    %   Reads melanoma data from file into an SCDataset class
    %   A known flaw with this dataset is that it contains double rows for the
    %   genes MARCH1 and MARCH2, a recommendation is to discard the
    %   duplicates, for example by ds = ds.geneSubset(unique(ds.genes));
    %   Publication: 'Dissecting the multicellular ecosystem of
    %   metastatic melanoma by single-cell RNA-seq', GSE72056
    %   Data is in TPM; this is not a UMI dataset.
    %
    % Johan Gustafsson, 2019-05-24
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsMel.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading melanoma...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/mel.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsMel.import('../../ImportableData/mel1/GSE72056_melanoma_single_cell_revised_v2.txt');
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
            %   path       Path to the 10x files. No slash at the end.
            %
            % Usage: ds = DsMel.import('../../ImportableData/mel1/GSE72056_melanoma_single_cell_revised_v2.txt');
            %
            
            ds = SCDataset;
            ds.name = 'mel';
            L = importdata(path,'\t');
            ds.data = L.data;
            [m,n] = size(ds.data);
            ds.cellIds = L.textdata(1, 2:(n+1));
            ds.genes = L.textdata(5:(m+1), 1);
            
            ds.sampleIds = ds.data(1, :);
            ds.sampleIds = sprintfc('%d',ds.sampleIds);
            malignant = ds.data(2, :);
            nonMalType = ds.data(3, :);
            ds.data = ds.data(4:m, :); %skip classification and tumor id in the data
            
            temp = malignant*10 + nonMalType; % create one number from 2 variables
            ds.cellType = arrayfun(@DsMel.transType, temp);
            
            % Transform to TPM according to (pow(2.0, dVal) - 1.0) * 10.0. Run TPM
            % afterwards just in case
            ds.data = (2.^ds.data-1)*10;
            ds = TPM(ds);
            ds = ds.fillEmpties();
        end
        
        function res = transType(x)
            if x <= 10
                res = Celltype.Unknown;
            elseif x >= 20
                res = Celltype.Malignant;
            else
                switch x
                    case 11
                        res = Celltype.TCell;
                    case 12
                        res = Celltype.BCell;
                    case 13
                        res = Celltype.Macrophage;
                    case 14
                        res = Celltype.Endothelial;
                    case 15
                        res = Celltype.Fibroblast;
                    case 16
                        res = Celltype.NKCell;
                    otherwise
                        res = Celltype.Unknown;
                end
            end
        end
    end
end

