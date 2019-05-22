classdef SCDep
    % SCDep
    %   Wrapper class for dataset access. The datasets are cached in
    %   persistant variables and also on disc, and only loaded when needed.
    %   To clear the persistant variables to get the memory back
    %   call "clear SCDep" from outside!
    %   It doesn't work to put that in a function, since the class is then
    %   used in the call and thereby locked
    %
    % Johan Gustafsson, 2019-05-21
    %
    
    properties
    end
    methods(Static)
        function init()
            % init
            %   Initialization code run once. Should be called in every
            %   function in this class.
            %
            % Usage: SCDep.init();
            %
            
            %set paths etc only once
            persistent inited
            if isempty(inited)
                inited = true;
                
                %put all .mat files in a temp folder, make sure it is created
                prevDir = SCDep.setPathToSource();
                if(~exist('../../TempData','dir'))
                    mkdir('../../TempData');
                end
                SCDep.restoreDir(prevDir);
            end
        end
        
        function prevDir = setPathToSource()
            % setPathToSource
            %   Sets the current directory to that of this source file. The
            %   old current directory is returned.
            %
            % Usage: oldDir = SCDep.setPathToSource();
            %
            codeDir = fileparts(which(mfilename));
            prevDir = pwd();
            cd(codeDir);
        end
        
        function restoreDir(prevDir)
            % restoreDir
            %   Restores the current directory after a call to
            %   setPathToSource(). 
            % Input:
            %   prevDir         Should be the return value from
            %                   setPathToSource
            % Usage: SCDep.restoreDir(oldDir);
            %
            cd(prevDir);
        end
        
        %Median from all gtex tissue samples
        function ret = gtex_mdn()
            SCDep.init();
            persistent v;
            if isempty(v)
                disp('reading GTEx median...');
                prevDir = SCDep.setPathToSource();
                filename = '../../TempData/GTExMed.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = ReadGTExMedian('../../ImportableData/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct');
                    save(filename, 'v');
                else
                    a = load(filename);
                    v = a.v;
                end
                SCDep.restoreDir(prevDir);
            end
            ret = v;
        end
        
        %ovarian ascites
        function ret = scd_ovasc()
            SCDep.init();
            persistent v;
            if isempty(v)
                disp('reading ovarian ascites data ...');
                prevDir = SCDep.setPathToSource();
                filename = '../../TempData/ovasc.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = ReadOvarianAscites('../../ImportableData/ovarian_ascites_data.txt', '../../ImportableData/ovarian_ascites_patids.txt', '../../ImportableData/ovarian_ascites_genes.txt', '../../ImportableData/classified_data_donor_abc_12k_merged.mat');
                    save(filename, 'v');
                else
                    a = load(filename);
                    v = a.v;
                end
                SCDep.restoreDir(prevDir);
            end
            ret = v;
        end
        
        %lung cancer 50000
        function [scd_lc_tumor,scd_lc_healthy] = scd_lc()
            SCDep.init();
            persistent v;
            persistent w;
            if isempty(v)
                disp('reading lung cancer data ...');
                prevDir = SCDep.setPathToSource();
                filename = '../../TempData/lc_tumor.mat';
                filenameh = '../../TempData/lc_healthy.mat';
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    [v,w] = ReadLung50000('../../ImportableData/LungCancer50000/AllCells_chunk1.txt', '../../ImportableData/LungCancer50000/AllCells_chunk2.txt', '../../ImportableData/LungCancer50000/AllCells_chunk3.txt', '../../ImportableData/LungCancer50000/MetaData.txt');
                    %make them sparse, both to save memory and to make it
                    %possible to save them
                    v.data = sparse(v.data);
                    w.data = sparse(w.data);
                    save(filename, 'v');
                    save(filenameh, 'w');
                else
                    a = load(filename);
                    v = a.v;
                    b = load(filenameh);
                    w = b.w;
                end
                SCDep.restoreDir(prevDir);
            end
            scd_lc_tumor = v;
            scd_lc_healthy = w;
        end
        
        %pbmc68000
        function ret = scd_pbmc68000()
            SCDep.init();
            persistent v;
            if isempty(v)
                disp('reading pbmc68000 data ...');
                filename = '../../TempData/pbmc68000.mat';
                prevDir = SCDep.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = ReadPBMC68000('../../ImportableData/PBMC68000PatAFresh/filtered_matrices_mex/hg19');
                    save(filename, 'v');
                else
                    a = load(filename);
                    v = a.v;
                end
                SCDep.restoreDir(prevDir);
            end
            ret = v;
        end
        
        %pbmcb10000
        function ret = scd_pbmcb10000()
            SCDep.init();
            persistent v;
            if isempty(v)
                disp('reading pbmcb10000 data ...');
                filename = '../../TempData/pbmcb10000.mat';
                prevDir = SCDep.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = ReadPBMCB10000('../../ImportableData/PBMC10000BCells/filtered_matrices_mex/hg19');
                    save(filename, 'v');
                else
                    a = load(filename);
                    v = a.v;
                end
                SCDep.restoreDir(prevDir);
            end
            ret = v;
        end
        
        %pbmctcd4mem10000
        function ret = scd_pbmctcd4mem10000()
            SCDep.init();
            persistent v;
            if isempty(v)
                disp('reading pbmccd4mem10000 data ...');
                filename = '../../TempData/pbmctcd4mem10000.mat';
                prevDir = SCDep.setPathToSource();
                if(~exist(filename,'file'))
                    disp('No .mat file found, importing data');
                    v = ReadPBMCTCD4Mem10000('../../ImportableData/PBMCCD4TCellsMemory/filtered_matrices_mex/hg19');
                    save(filename, 'v');
                else
                    a = load(filename);
                    v = a.v;
                end
                SCDep.restoreDir(prevDir);
            end
            ret = v;
        end
        
    end
    
    
end

