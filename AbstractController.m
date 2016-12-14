classdef (Abstract) AbstractController < handle
    % This is the generic controller
    
    properties
        controlled_link_ids % [Ix1] ids of controlled ramps
        control_sequence @ ControlSequence
    end
    
    methods
        function this = AbstractController(model)
            this.controlled_link_ids = Utils.column(model.controlled_link_ids);
            this.control_sequence = ControlSequence(this, this.controlled_link_ids);
        end
    end
    
    methods (Abstract)
        
        % returns a singleton ControlSequence
        [u] = get_control(this,state,time);
        
        % returns a ControlSequence
        [u] = compute_control_sequence(this,initial_state,predicted_demands);
        
    end
    
end

