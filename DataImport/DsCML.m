classdef DsCML
    % DsCML
    %   Reads CML data from file into an SCDataset class
    %   The cells are CD38- CD34+ stem cells, checked for BCR-ABL
    %   BCR-ABL-positive cells are of the cell type malignant in the
    %   dataset, even though also stem cells.
    %   Pooled normal stem cells, will have 'healthy bulk' in the 
    %   extraCellInfo field
    %   Cell line cells are marked with 'cell line' in the extraCellInfo
    %   field.
    %   There is additional info about CML progression etc. that we
    %   currently do not import.
    %   Publication: 'Single-cell transcriptomics uncovers distinct 
    %   molecular signatures of stem cells in chronic myeloid leukemia',
    %   GSE76312.
    %   Data is in TPM; this is not a UMI dataset.
    %
    % Johan Gustafsson, 2019-05-25
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsCML.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading CML...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/cml.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsCML.import('../../ImportableData/GSE76312_CML/Giustacchini_Thongjuea_et.al_Nat.Med.RPKM.txt', '../../ImportableData/GSE76312_CML/GSE76312_Giustacchini-Thongjuea_et.al_Cell_Annotation.txt');
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
        function  ds = import(path_, pathMetadata)
            % import
            %   Imports the data
            % Input:
            %   path_           Path to the data
            %   pathMetadata    Path to the metadata file
            %
            % Usage: ds =
            % DsCML.import('../../ImportableData/GSE76312_CML/Giustacchini_Thongjuea_et.al_Nat.Med.RPKM.txt',
            %              '../../ImportableData/GSE76312_CML/GSE76312_Giustacchini-Thongjuea_et.al_Cell_Annotation.txt');
            %
            
            ds = SCDataset;
            ds.name = 'cml';
            L = importdata(path_,'\t');
            ds.data = TPM(L.data);
            [m,n] = size(ds.data);
            ds.cellIds = strsplit(L.textdata{1,1});
            ds.genes = L.textdata(2:end,1);
            
            %get metadata
            t = readtable(pathMetadata, 'ReadVariableNames',true, 'ReadRowNames', false, 'Delimiter', '\t');
            
            %I checked that the cell ids comes in the same order, so no
            %need to do any mapping between ids
            
            ds = ds.fillEmpties();
            
            %The BCR_ABL_status field can have the following values:
            %'negative' - cellType: stem cell
            %'negative_low_gapdh' - cellType: Unknown
            %'negative_undetected_gapdh' - cellType: Unknown
            %'no_primers' - Cell line, will be removed in the dataset
            %'normal' - These are pooled normal stem cells, so not single cells. 
            %           Will have 'healthy bulk' in the extraCellInfo field
            %'positive' - cellType: malignant
            %'primers' - Cell line, will be removed in the dataset
            bcrabl = t.BCR_ABL_status.';
            isMalignant = strcmp(bcrabl,'positive');
            isUnknown = strcmp(bcrabl,'negative_undetected_gapdh') | strcmp(bcrabl,'negative_undetected_gapdh');
            isCellLine = strcmp(bcrabl,'primers') | strcmp(bcrabl,'no_primers');
            isBulk = strcmp(bcrabl,'normal');
            ds.cellType(1,:) = Celltype.HematopeticStemOrProgenitor;
            ds.cellType(1,isMalignant) = Celltype.Malignant;
            ds.cellType(1,isUnknown) = Celltype.Unknown;
            ds.extraCellInfo(1,isBulk) = {'healthy bulk'};
            ds.extraCellInfo(1,isCellLine) = {'cell line'};
            
            ds.sampleIds = t.Patient_id.';
            
        end
    end
end

