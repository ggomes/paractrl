classdef LinearizedMPCController < PredictiveController
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function this = LinearizedMPCController(num_base_controllers)
            this@PredictiveController(num_base_controllers);
        end
           
        function [u]=compute_control_sequence(this,warm_start,initial_state,predicted_demands)
            u = ControlSequence(this,rand(1,length(warm_start.control_sequence)));
        end
        
    end
    
end

