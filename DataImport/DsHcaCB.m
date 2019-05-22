classdef DsHcaCB
    % DsHcaCB
    %   Reads HCA cord blood from file into an SCDataset.
    %   The dataset covers roughly 354000 cells.
    %   Classification was received from the authors (Bo Li)
    %
    % Johan Gustafsson, 2019-05-22
    %
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsHcaCB.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading hca cord blood...');
                filename = '../../TempData/hca_cb.mat';
                filename2 = '../../TempData/hca_cb_2.mat';
                filename3 = '../../TempData/hca_cb_3.mat';
                prevDir = DsHelper.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsHcaCB.import('../../ImportableData/ica_cord_blood_h5.h5', '../../ImportableData/cord_blood_cell_type_processed.txt');
                    %cut up the data in several blocks and save to separate files
                    %since it seems there is a limitation to object size
                    data2 = v.data(:,100001:200000);
                    data3 = v.data(:,200001:end);
                    v.data = v.data(:,1:100000);
                    
                    save(filename, 'v');
                    save(filename2, 'data2');
                    save(filename3, 'data3');
                else
                    a = load(filename);
                    v = a.v;
                    b = load(filename2);
                    c = load(filename3);
                    data2 = b.data2;
                    data3 = c.data3;
                end
                %put the data matrix back together
                v.data = [v.data data2 data3];
                
                DsHelper.restoreDir(prevDir);
            end
            ret = v;
        end
    end
    
    methods(Static, Access = private)
        
        function ds = import(path, classPath)
            % import
            %   Imports the data.
            % Input:
            %   path                Path to the .h5 file
            %   classPath           Path to the cell type info file
            %
            % Usage: ds = DsHcaCB.import('../../ImportableData/ica_cord_blood_h5.h5', '../../ImportableData/cord_blood_cell_type_processed.txt');
            %
            
            %path = '../../ImportableData/ica_cord_blood_h5.h5';
            %path = 'C:/Work/MatlabCode/components/SCLib/ImportableData/ica_cord_blood_h5.h5';
            %classPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/HCA_cb_ct.txt';
            data = h5read(path,'/GRCh38/data');
            %genes = h5read(path,'/GRCh38/genes');%this is ensembl gene ids
            geneNames = h5read(path,'/GRCh38/gene_names');
            indices = h5read(path,'/GRCh38/indices');
            indptr = h5read(path,'/GRCh38/indptr');
            shape = h5read(path,'/GRCh38/shape');
            barcodes = h5read(path,'/GRCh38/barcodes');
            %create a vector with cell index for each value
            numCells = shape(2);
            numDataPoints = size(data, 1);
            dataPointCellIndex = zeros(numDataPoints, 1);
            lastPoint = 0;%zero indexed
            for i = 1:numCells
                %everything is zero indexed in the data file which makes this tricky...
                nextPoint = indptr(i+1);
                dataPointCellIndex(lastPoint+1:nextPoint) = i;
                lastPoint = nextPoint;
            end
            
            %remove a bunch of char(0) at the end:
            barcodes = regexprep(barcodes,char(0),'');
            
            
            ds = SCDataset;
            ds.name = 'hca_cord';
            matches = regexp(geneNames, '[a-zA-Z0-9_\-\.]*', 'match');
            ds.genes = cellfun(@(x) x{1}, matches,'UniformOutput',false);
            ds.cellIds = barcodes.';
            %extract sample ids from cell ids
            %typical string: MantonCB1_HiSeq_1-AAAGCAACACTTGGAT-1
            ds.sampleIds = cellfun(@(x) x(7:9), ds.cellIds, 'UniformOutput', false);
            
            
            
            numGenes = size(ds.genes,1);
            
            dataPointGeneIndex = zeros(numDataPoints, 1);
            dataPointGeneIndex(:,1) = indices(:,1) + 1;
            dataVector = zeros(numDataPoints, 1);
            dataVector(:,1) = data(:,1);
            ds.data = sparse(dataPointGeneIndex, dataPointCellIndex, dataVector, double(numGenes), double(numCells));
            
            ds = ds.fillEmpties();
            
            %read classification from authors (it has been purified with some regexp replaces in Visual Studio, see file under scripts)
            f = readtable(classPath, 'ReadRowNames', false, 'Delimiter', '\t');
            c = table2cell(f);
            %transform from string type to SCLib cell type enum
            ct = cellfun(@DsHcaCB.String2CellTypeId, c(:,2));
            
            [~, ia, ib] = intersect(ds.cellIds, c(:, 1));
            ds.cellType(1,ia) = ct(ib,1).';
            
            %remove the cells with cell type 'unknown' - these are low quality cells
            %that did not have a classification in the file
            ds = ds.cellSubset(ds.cellType ~= Celltype.Unknown);
            
            
        end
        
        function ret = String2CellTypeId(str)
            str = strtrim(str);
            if strcmp(str,'CD4+ Naive T cells')
                ret = Celltype.TCellCD4Pos;%Skipped the naive part for now
            elseif strcmp(str,'CD8+ Naive T cells')
                ret = Celltype.TCellCD8Pos;%Skipped the naive part for now
            elseif strcmp(str,'T cells')
                ret = Celltype.TCell;
            elseif strcmp(str,'Naive T cells')
                ret = Celltype.TCell;%Skipped the naive part for now
            elseif strcmp(str,'B cells')
                ret = Celltype.BCell;
            elseif strcmp(str,'CD14+ Monocytes')
                ret = Celltype.Monocyte;%ignore CD14+ for now
                %elseif strcmp(str,'FCGR3A+ Monocytes')
                %    ret = Celltype.Monocyte;
            elseif strcmp(str,'NK cells')
                ret = Celltype.NKCell;
            elseif strcmp(str,'cDCs')
                ret = Celltype.Dendritic;
                %elseif strcmp(str,'Megakaryocytes')
                %    ret = Celltype.Megakaryocyte;
            elseif strcmp(str,'HSCs')
                ret = Celltype.HematopeticStemOrProgenitor;
            elseif strcmp(str,'Erythrocytes')
                ret = Celltype.Erythrocyte;
            else
                disp(strcat('unknown type: ',str));
                ret = Celltype.Unknown;
            end
        end
    end
end