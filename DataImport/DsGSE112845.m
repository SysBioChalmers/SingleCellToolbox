classdef DsGSE112845
    % DsBC2
    %   Reads one of the GSE112845 datasets from file into an SCDataset. Note
    %   that there are several datasets to choose from.
    % Input:
    %   path                Path to the 10x data folder
    %   classificationPath  Path to the cell type info file, can be omitted
    %
    % Usage: ds = ReadGSE112845('../../ImportableData/GSE112845/DTM-X_PBMC_live', '../../ImportableData/GSE112845/DTM-X_PBMC_live_ct.txt');
    %
    % Johan Gustafsson, 2019-05-20
    %
    methods(Static)
        function [xPat,yPat,cd8] = get()
            % get
            %   Gets the datasets. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: [xPat,yPat,cd8] = DsGSE112845.get();
            DsHelper.init();
            persistent v;
            persistent w;
            persistent x;
            if isempty(v)
                disp('reading GSE112845 data ...');
                filename = '../../TempData/gse112845_1.mat';
                filename2 = '../../TempData/gse112845_2.mat';
                filename3 = '../../TempData/gse112845_3.mat';
                prevDir = DsHelper.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsGSE112845.import('../../ImportableData/GSE112845/DTM-X_PBMC_live', '../../ImportableData/GSE112845/DTM-X_PBMC_live_ct.txt');
                    w = DsGSE112845.import('../../ImportableData/GSE112845/DTM-Y_PBMC_methanol_3w', '../../ImportableData/GSE112845/DTM-Y_PBMC_methanol_3w_ct.txt');
                    x = DsGSE112845.import('../../ImportableData/GSE112845/CD8_live', '');
                    v.name = 'GSE112845 DTM-X_PBMC_live';
                    w.name = 'GSE112845 DTM-Y_PBMC_methanol_3w';
                    x.name = 'GSE112845 CD8_live';
                    save(filename, 'v');
                    save(filename2, 'w');
                    save(filename3, 'x');
                else
                    a = load(filename);
                    v = a.v;
                    b = load(filename2);
                    w = b.w;
                    c = load(filename3);
                    x = c.x;
                end
                DsHelper.restoreDir(prevDir);
            end
            xPat = v;
            yPat = w;
            cd8 = x;
        end
    end
    
    methods(Static, Access = private)
        
        
        function ds = import(path, classificationPath)
            % import
            %   Reads one of the GSE112845 datasets from file into an SCDataset. Note
            %   that there are several datasets to choose from.
            % Input:
            %   path                Path to the 10x data folder
            %   classificationPath  Path to the cell type info file, can be omitted
            %
            % Usage: ds = ReadGSE112845('../../ImportableData/GSE112845/DTM-X_PBMC_live', '../../ImportableData/GSE112845/DTM-X_PBMC_live_ct.txt');
            %
            % Johan Gustafsson, 2019-05-20
            %
            
            %path = 'C:/Work/MatlabCode/components/SCLib/ImportableData/GSE112845/DTM-X_PBMC_live';
            %classificationPath = 'C:/Work/MatlabCode/components/SCLib/ImportableData/GSE112845/DTM-X_PBMC_live_ct.txt';
            
            ds = Read10xMatrix(path);
            ds.name = 'GSE112845';
            
            if ~isempty(classificationPath)
                %read ds.cellType file from the authors
                f = readtable(classificationPath, 'ReadRowNames', false, 'Delimiter', '\t');
                c = table2cell(f);
                %transform from string type to SCLib cell type enum
                ct = cellfun(@DsGSE112845.String2CellTypeId, c(:,2));
                
                [~, ia, ib] = intersect(ds.cellIds, c(:, 1));
                ds.cellType(1,ia) = ct(ib,1).';
                
                %remove the cells with cell type 'unknown' - these are low quality cells
                %that did not have a classification in the file
                ds = ds.cellSubset(ds.cellType ~= Celltype.Unknown);
            else
                %assume it is the CD8+ dataset
                ds.cellType(1,:) = Celltype.TCellCD8Pos;
                ds.cellType = ds.cellType;
            end
        end
        
        function ret = String2CellTypeId(str)
            if strcmp(str,'CD4 T cells')
                ret = Celltype.TCell;%subclass of t cells not reliable
            elseif strcmp(str,'CD8 T cells')
                ret = Celltype.TCell;
            elseif strcmp(str,'B cells')
                ret = Celltype.BCell;
            elseif strcmp(str,'CD14+ Monocytes')
                ret = Celltype.Monocyte;%ignore subset for now
            elseif strcmp(str,'FCGR3A+ Monocytes')
                ret = Celltype.Monocyte;
            elseif strcmp(str,'NK cells')
                ret = Celltype.NKCell;
            elseif strcmp(str,'Dendritic cells')
                ret = Celltype.Dendritic;
            elseif strcmp(str,'Megakaryocytes')
                ret = Celltype.Megakaryocyte;
            else
                disp(strcat('unknown type: ',str));
                ret = Celltype.Unknown;
            end
        end
        
    end
end