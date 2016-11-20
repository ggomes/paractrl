classdef DemandProfile < handle
    
    properties
        time        % 1 x K (seconds after midnight)
        link_ids    % I x 1 (source links)
        demands     % I x K (veh / time step, e.g. veh every 4 seconds)
    end
    
    methods
        
        function this = DemandProfile(link_ids)
            this.link_ids = link_ids;
            this.time = [];
            this.demands = [];
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
        
        function []=plot(this)
            figure
            plot(this.time,this.demands)
            grid
            xlabel('time [sec]')
            ylabel('flow [vph]')
        end
    end
    
end

