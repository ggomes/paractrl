classdef State
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        link_ids    % I x 1
        link_dty    % I x 1 [veh / link]
    end
    
    methods
        
        function this = State(link_ids)
            this.link_ids = link_ids;
            this.link_dty = [];
        end
        
    end
    
end

