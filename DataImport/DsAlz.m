classdef DsAlz
    % DsAlz
    %   Reads the Alzheimer's cells from file into an SCDataset.
    %   The dataset covers roughly 70000 cells from 48 patients.
    %   
    %
    % Johan Gustafsson, 2019-06-02
    %
    
    %For now, the dataset returned sets cell type to a string instead of a
    %number, since it is just easier; it seems like the cell types are very
    %specialized.
    
    
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsAlz.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading Alzheimer data ...');
                filename = '../../TempData/alz.mat';
                prevDir = DsHelper.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsAlz.import('../../ImportableData/Alzheimers');
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
    
    methods(Static, Access = public)
        function ds = import(directoryPath)
            % import
            %   Imports the data.
            % Input:
            %   directoryPath       Path to the 10x files. No slash at the end.
            %
            % Usage: ds = import('../../ImportableData/Alzheimer');
            %
            
            
            %classificationPath = strcat(directoryPath,'/68k_pbmc_barcodes_annotation.tsv');
            ds = Read10xMatrix(directoryPath, ...
                               'filtered_count_matrix.mtx', ...
                               'filtered_gene_row_names.txt', ...
                               'filtered_column_metadata.txt', ...
                               true, 1);
            ds.name = 'Alz';
            
            %read the metadata file once more to get the cell types and 
            f = readtable(strcat(directoryPath,'/filtered_column_metadata.txt'), 'ReadRowNames', false, 'Delimiter', '\t', 'FileType', 'text');
            
            c = table2cell(f);
            ds.sampleIds = f.projid.';
            %ct = cellfun(@DsAlz.String2CellTypeId, f.broad_cell_type.');
            ds.extraCellInfo = f.Subcluster.';
            
            [~, ia, ib] = intersect(ds.cellIds.', c(:, 1));
            %ds.cellType(1,ia) = ct(ib,1).';
            %make cell type string for now to simplify things
            ds.cellType = ds.extraCellInfo;%just an empty vector to allow for next line
            ds.cellType(1,ia) = f.broad_cell_type(ib,1).';
            ds.extraCellInfo(1,ia) = f.Subcluster(ib,1).';
            
        end
        %not used for now
        %{
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
        %}
    end
end
