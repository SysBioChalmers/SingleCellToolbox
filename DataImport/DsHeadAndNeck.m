classdef DsHeadAndNeck
    % DsHeadAndNeck
    %   Reads head and neck carcinoma data from file into an SCDataset class
    %   Publication: 'Single-Cell Transcriptomic Analysis of Primary and 
    %   Metastatic Tumor Ecosystems in Head and Neck Cancer', GSE103322
    %   Data is in TPM; this is not a UMI dataset.
    %   We had to remove row 6 (the classifications, as strings...) from the 
    %   file and put it in a separate file to be able to read the data in 
    %   an easy way...
    %
    % Johan Gustafsson, 2020-01-21
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsHeadAndNeck.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading HeadAndNeck...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/HeadAndNeck.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsHeadAndNeck.import('../../ImportableData/HANC/HNSCC_all_datamod.txt', '../../ImportableData/HANC/HNSCC_Separate_Classification.txt');
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
        function  ds = import(path, pathClass)
            % import
            %   Imports the data
            % Input:
            %   path       Path to the data file
            %   pathClass  Path to cell classifications file
            %
            % Usage: ds = DsHeadAndNeck.import('../../ImportableData/HNSCC_all_datamod.txt', '../../ImportableData/HNSCC_Separate_Classification.txt');
            %
            
            ds = SCDataset;
            ds.name = 'hanc';
            L = importdata(path,'\t');
            [m,n] = size(L.data);
            [r,s] = size(L.textdata);
            ds.data = L.data(5:m,:);
            ds.cellIds = L.textdata(1, 2:(s));
            ds.genes = L.textdata(6:r, 1);
            %get rid of the quotes in the gene names
            matches = regexp(ds.genes, '[a-zA-Z0-9_\-]*', 'match');
            ds.genes = cellfun(@(x) x{1}, matches,'UniformOutput',false);

            %for classification, assume that it is enough to look at the "classified as
            %cancer cell" and "non-cancer cell type"
            %read the classification from the separate file

            fileID = fopen(pathClass,'r');
            line = fgetl(fileID);
            types = strsplit(line,'\t');
            fclose(fileID);
            [t, u] = size(types);
            types = types(2:u);
            malignant = L.data(3, :);
            nonMalType = arrayfun(@DsHeadAndNeck.String2CellTypeId, types,'UniformOutput',false);
            nonMalType = cell2mat(nonMalType).';
            cc = logical(malignant);
            ds.cellType = zeros(1, n);
            ds.cellType(1, cc) = 2;%tumor
            ds.cellType(1, ~cc) = nonMalType(~cc);

            %read tumor ids
            matches = regexp(ds.cellIds, '^[a-zA-Z0-9]+', 'match');
            ds.sampleIds = cellfun(@(x) x{1}, matches,'UniformOutput',false);


            % Transform to TPM according to (pow(2.0, dVal) - 1.0) * 10.0, and run TPM
            % just to make sure
            ds.data = (2.^ds.data-1)*10;
            ds = TPM(ds);
        end
        
        function ret = String2CellTypeId(str)
            if strcmp(str,'Fibroblast')
                ret = Celltype.Fibroblast;
            elseif strcmp(str,'-Fibroblast')
                ret = Celltype.Fibroblast;
            elseif strcmp(str,'B cell')
                ret = Celltype.BCell;
            elseif strcmp(str,'myocyte')
                ret = Celltype.Myocyte;
            elseif strcmp(str,'T cell')
                ret = Celltype.TCell;
            elseif strcmp(str,'Dendritic')
                ret = Celltype.Dendritic;
            elseif strcmp(str,'Mast')
                ret = Celltype.Mast;
            elseif strcmp(str,'Macrophage')
                ret = Celltype.Macrophage;
            elseif strcmp(str,'Endothelial')
                ret = Celltype.Endothelial;
            elseif strcmp(str,'0')
                ret = Celltype.Unknown;
            else
                disp('ReadHeadAndNeck: Failed to convert celltype string');
                disp(str);
                ret = Celltype.Unknown;
            end
        end
    end
end

