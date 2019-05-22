classdef DsPbmc68k
    % DsPbmc68k
    %   Reads the PBMC68k cells from file into an SCDataset.
    %   The dataset covers roughly 68000 cells from blood.
    %   The file 68k_pbmc_barcodes_annotation.tsv is expected in directoryPath
    %
    % Johan Gustafsson, 2019-05-22
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsOvAsc.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading pbmc68000 data ...');
                filename = '../../TempData/pbmc68000.mat';
                prevDir = DsHelper.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsPbmc68k.import('../../ImportableData/PBMC68000PatAFresh/filtered_matrices_mex/hg19');
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
        function ds = import(directoryPath)
            % import
            %   Imports the data.
            % Input:
            %   directoryPath       Path to the 10x files. No slash at the end.
            %
            % Usage: ds = import('../../ImportableData/PBMC68000PatAFresh/filtered_matrices_mex/hg19');
            %
            
            
            %directoryPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/PBMC68000PatAFresh/filtered_matrices_mex/hg19';
            classificationPath = strcat(directoryPath,'/68k_pbmc_barcodes_annotation.tsv');
            
            ds = Read10xMatrix(directoryPath);
            ds.name = 'pbmc 68000';
            
            f = readtable(classificationPath, 'ReadRowNames', false, 'Delimiter', '\t', 'FileType', 'text');
            c = table2cell(f);
            ct = cellfun(@DsPbmc68k.String2CellTypeId, c(:,4));
            
            [~, ia, ib] = intersect(ds.cellIds.', c(:, 3));
            ds.cellType(1,ia) = ct(ib,1).';
            
        end
        
        function ret = String2CellTypeId(str)
            %we ignore all detailed subsets for now
            if strcmp(str,'CD8+ Cytotoxic T')
                ret = Celltype.TCellCD8Pos;
            elseif strcmp(str,'CD8+/CD45RA+ Naive Cytotoxic')
                ret = Celltype.TCellCD8Pos;
            elseif strcmp(str,'CD4+/CD45RO+ Memory')
                ret = Celltype.TCellCD4Pos;
            elseif strcmp(str,'CD4+/CD45RA+/CD25- Naive T')
                ret = Celltype.TCellCD4Pos;
            elseif strcmp(str,'CD4+ T Helper2')
                ret = Celltype.TCellCD4Pos;
            elseif strcmp(str,'CD4+/CD25 T Reg')
                ret = Celltype.TCellReg;
            elseif strcmp(str,'CD19+ B')
                ret = Celltype.BCell;
            elseif strcmp(str,'CD14+ Monocyte')
                ret = Celltype.Monocyte;
            elseif strcmp(str,'CD56+ NK')
                ret = Celltype.NKCell;
            elseif strcmp(str,'Dendritic')
                ret = Celltype.Dendritic;
            elseif strcmp(str,'CD34+')
                ret = Celltype.HematopeticStemOrProgenitor;
                %elseif strcmp(str,'Megakaryocytes')
                %    ret = Celltype.Megakaryocyte;
            else
                disp(strcat('unknown type: ',str));
                ret = Celltype.Unknown;
            end
        end
    end
end
