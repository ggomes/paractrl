classdef Utils
    
    methods(Static)
        
        function [c]=column(x)
            c = reshape(x,numel(x),1);
        end
        
        function [r]=row(x)
            r = reshape(x,1,numel(x));
        end
    end
    
end

