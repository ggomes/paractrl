classdef BaseBlock < handle
    %UNTITLED14 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        controllers 
    end
    
    methods
        
        function this = BaseBlock(controllers)
            this.controllers = controllers;
        end
        
        function [ctrl_seq_set] = compute_control_sequence(this,initial_state,predicted_demands)
            
            ctrl_seq_set = [];
            for i=1:length(this.controllers)
                this.controllers{i}.compute_control_sequence(initial_state,predicted_demands);
                ctrl_seq_set = [ctrl_seq_set this.controllers{i}.control_sequence];
            end
            
        end
        
        function [n]=num_controllers(this)
            n = length(this.controllers);
        end
        
    end
    
end

