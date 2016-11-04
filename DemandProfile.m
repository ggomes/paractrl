classdef DemandProfile
    
    properties
        time 
        values_map@containers.Map
    end
    
    methods
        
        function this = DemandProfile(time)
            this.time = time;
            this.values_map = containers.Map('KeyType','int32','ValueType','any');
        end
        
        function this = add_values(this,link_id,values)
            if length(values)~=length(this.time)
                error('length(values)~=length(this.time)')
            end
            this.values_map(link_id) = values;
        end
        
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

