classdef StateTrajectory < handle
    
    properties
        time        % 1 x K [seconds after midnight]
        link_ids    % I x 1
        link_dty    % I x K [veh / link]
    end
    
    methods
        
        function this = StateTrajectory(link_ids)
            this.link_ids = link_ids;
            I = length(link_ids);
            this.link_dty = [];
            this.time = [];
        end
        
        % NOTE DOCUMENT THIS!
        function this = add_values(this,t,values)
            this.time = [this.time t];
            this.link_dty = [this.link_dty values];
        end
        
    end
    
end

