classdef SimpleController < AbstractController
    
    properties
        rate        % [units?]
    end
    
    methods
        
        function this = SimpleController(model,rate)
            this@AbstractController(model);
            this.rate = rate;
        end
        
        function [] = compute_control_sequence(this,initial_state,predicted_demands)
            this.control_sequence = ControlSequence(this,this.source_link_ids);
            this.control_sequence.time = predicted_demands.time;
            this.control_sequence.control_sequence = this.rate*ones(1,numel(predicted_demands.time));
        end
        
        function [u] = get_control(this,state,time)
            u = this.rate;
        end
        
    end
    
end

