classdef DemandProfile < handle
    
    properties
        time        % 1 x K (seconds)
        link_ids    % I x 1 (source links)
        demands     % I x K (veh / time step)
    end
    
    methods
        
        function this = DemandProfile(link_ids)
            this.link_ids = link_ids;
            this.time = [];
            this.demands = [];
        end
        
        % NOTE DOCUMENT THIS
        function this = add_values(this,t,d)
            this.time = [this.time t];
            this.demands = [this.demands d];
        end
        
        % NOTE DOCUMENT THIS
        function this = perturb(this,delta)
            warning('implement this')
        end
        
        function [x] = get_initial_time(this)
            x = this.time(1);
        end
        
        function [x] = get_final_time(this)
            x = this.time(end);
        end
        
    end
    
end

