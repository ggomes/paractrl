classdef (Abstract) AbstractController < handle
    % This is the generic controller
    
    properties
        control_sequence % this is an array of control_sequences
                         % computed by compute_control_sequence
    end
    
    methods
        function this = AbstractController(num_control_sequence)
            this.control_sequence = ControlSequence.empty(num_control_sequence,0);
        end
    end
    
    methods (Abstract)
        [u] = get_control(this,state,time);
        
        % returns a control sequence
        [u] = compute_control_sequence(this,initial_state,predicted_demands);
    end
    
end

