classdef DsHelper
    % DsHelper
    %   Helper class for all dataset import classes. In general these classes 
    %   support caching of data, both at
    %   persistant variables level but also on disc, and only loaded when needed.
    %   To clear the persistant variables to get the memory back
    %   call "clear <class name>" from outside!
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
            % Usage: DsBase.init();
            %
            
            %set paths etc only once
            persistent inited
            if isempty(inited)
                inited = true;
                
                %put all .mat files in a temp folder, make sure it is created
                prevDir = DsHelper.setPathToSource();
                if(~exist('../../TempData','dir'))
                    mkdir('../../TempData');
                end
                DsHelper.restoreDir(prevDir);
            end
        end
        
        function prevDir = setPathToSource()
            % setPathToSource
            %   Sets the current directory to that of this source file. The
            %   old current directory is returned.
            %
            % Usage: oldDir = DsHelper.setPathToSource();
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
            % Usage: DsHelper.restoreDir(oldDir);
            %
            cd(prevDir);
        end
    end
    
    
end

