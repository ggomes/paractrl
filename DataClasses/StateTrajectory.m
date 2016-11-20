classdef StateTrajectory < handle
    
    properties
        time        % 1 x K [seconds after midnight]
        link_ids    % I x 1
        link_veh    % I x K [veh / link]
    end
    
    methods
        
        function this = StateTrajectory(link_ids)
            this.link_ids = link_ids;
            this.link_veh = [];
            this.time = [];
        end
        
        % NOTE DOCUMENT THIS!
        function this = add_values(this,t,state)
            this.time = [this.time t];
            this.link_veh = [this.link_veh state.link_veh];
        end
        
        function this = plot(this)
           
            figure
            h = pcolor(this.link_veh);
            set(h,'EdgeAlpha',0)
            
        end
        
    end
    
end

