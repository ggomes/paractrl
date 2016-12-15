classdef ControlSequence < handle
    
    properties
        controller @ AbstractController
        time       % 1 x K [in seconds after midnight]
        link_ids   % I x 1  (controlled ramps) 
        rate_vph   % I x K  [vph]
    end
    
    methods
        
        function this = ControlSequence(controller,link_ids)
            this.controller = controller;
            this.link_ids = Utils.column(link_ids);
            this.rate_vph = nan(numel(link_ids),0);
            this.time = [];
        end
        
        % NOTE DOCUMENT THIS!
        function this = add_values(this,time,rate_vph)
            this.time = [this.time time];
            this.rate_vph = [this.rate_vph rate_vph];            
        end
        
        function [x] = is_singleton(this)
            x = numel(this.time)==1;
        end
        
    end
    
end

