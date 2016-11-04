classdef ParametrizedMPCController < PredictiveController
    %UNTITLED10 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        warmup_parameters
    end
    
    methods
        
        function this = ParametrizedMPCController(num_base_controllers)
            this@PredictiveController(num_base_controllers);
        end
        
        function this = set_warm_start(this,u)
            this.set_warm_start@PredictiveController(u);
            this.warmup_parameters  = [1 2 3]; % translate_parameters(u);
        end
           
        function [u]=compute_control_sequence(this,warm_start,initial_state,predicted_demands)
            u = ControlSequence(this,rand(1,length(warm_start.control_sequence)));
        end
        
    end
    
end

