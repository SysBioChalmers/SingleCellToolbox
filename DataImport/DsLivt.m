classdef DsLivt
    % DsLivt
    %   Reads Liver cancer T cells from file into an SCDataset.
    %   The dataset covers roughly 4000 cells.
    %   No classification other than that they are T cells is available.
    %
    % Johan Gustafsson, 2019-05-20
    %
    
    methods(Static)
        function ret = get()
            % get
            %   Gets the dataset. This is quick except the first time it is
            %   called, since the data is cached at two levels; a
            %   persistant variable and in a .mat file.
            % Usage: ds = DsLivt.get();
            DsHelper.init();
            persistent v;
            if isempty(v)
                disp('reading liver cancer T cells...');
                prevDir = DsHelper.setPathToSource();
                filename = '../../TempData/livt.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = DsLivt.import('../../ImportableData/GSE98638_HCC.TCell.S5063.count.txt', '../../ImportableData/GSE98638_HCC.TCell.OkCellIds.txt');
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
        function  ds = import(pathData, pathCellFilter)
            % import
            %   Imports the data.
            % Input:
            %   path                Path to the data file
            %   pathCellFilter      Path to the cell quality file
            %
            % Usage: ds = DsLivt.import('../../ImportableData/GSE98638_HCC.TCell.S5063.count.txt', '../../ImportableData/GSE98638_HCC.TCell.OkCellIds.txt');
            %
            
            ds = SCDataset;
            ds.name = 'LCTCells';
            L = importdata(pathData,'\t');
            ds.data = L.data;
            [m,n] = size(ds.data);
            ds.cellIds = L.textdata(1, 3:end);
            ds.genes = L.textdata(2:end, 2);
            %filter on cell filter, i.e. remove the cells that did not pass quality
            %control in the paper:
            f = readtable(pathCellFilter, 'ReadVariableNames',false, 'ReadRowNames', false, 'Delimiter', '\t');
            okCells = table2cell(f);
            [ds.cellIds, ia, ib] = intersect(ds.cellIds, okCells);
            ds.data = ds.data(:, ia);
            
            %extract patient id, format is either 'xxxx-xx-yyyy' or 'xxxx-yyyy', and we
            %only want the y:s
            %a bit annoying, we have to extract two tokens and only keep the last
            temp = regexp(ds.cellIds, '\w+-(\w+-)?(\w+)', 'tokens');
            ds.sampleIds = cellfun(@(c) c{1}{2},temp, 'UniformOutput', false);
            ds.cellType = repmat(Celltype.TCell,size(ds.genes,1),size(ds.cellIds,2));
            
            ds = ds.fillEmpties();
        end
    end
end