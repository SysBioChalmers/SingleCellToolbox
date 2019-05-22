classdef ProgrBarContext
    % ProgrBarContext
    %   Support class for the ProgrBarContext class. Should normally be
    %   constructed via functions in the ProgrBar, not directly, even
    %   though it is possible if desired. The class represents a progress
    %   bar context, which can be empty (no parent), silent or with a parent.
    %
    % Johan Gustafsson, 2019-05-21
    %
    properties
        silent % if true no output will be made
        parent % parent progress bar
        fraction % fraction of the parent's progress bar
    end
    methods
        function obj = ProgrBarContext(silent_, parent_, fraction_)
            % ProgrBarContext
            %   Constructor for the ProgrBarContext class.
            % Input:
            %   silent_     (optional) True if the context should be
            %               silent. Defaults to false.
            %   parent_     (optional) The parent progress bar. Could
            %               either be [] or a ProgrBar object (since
            %               ProgrBar is a handle class, it will be a
            %               reference to a ProgrBar object). Defalts to [].
            %   fraction_   (optional) The fraction of the parents total 
            %               progress that is allocated to this context.
            %               Defaults to 0.
            % Usage: ctxt = ProgrBarContext(false, bar, 0.3);
            %
            if nargin < 1
                silent_ = false;
            end
            if nargin < 2
                parent_ = [];
            end
            if nargin < 3
                fraction_ = 0;
            end
            
            obj.silent = silent_;
            obj.parent = parent_;
            obj.fraction = fraction_;
        end
    end
end
    