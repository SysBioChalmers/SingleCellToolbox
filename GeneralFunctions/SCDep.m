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
        
        
    end
    
    
end

