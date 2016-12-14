classdef ControlSequence < handle
    
    properties
        controller @ AbstractController
        time                 % 1 x K [in seconds after midnight]
        link_ids             % I x 1  (controlled ramps) 
        control_sequence     % I x K  [veh/time step, e.g. vehicles every 4 seconds]
    end
    
    methods
        
        function this = ControlSequence(controller,link_ids)
            this.controller = controller;
            this.link_ids = Utils.column(link_ids);
            this.control_sequence = nan(numel(link_ids),0);
            this.time = [];
        end
        
        % NOTE DOCUMENT THIS!
        function this = add_values(this,t,c)
            this.time = [this.time t];
            this.control_sequence = [this.control_sequence c];            
        end
        
        function [x] = is_singleton(this)
            x = numel(this.time)==1;
        end
        
    end
    
end

