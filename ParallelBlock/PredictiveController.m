classdef PredictiveController < AbstractController
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods
        
        function this = PredictiveController(num_base_controllers)
            this@AbstractController(num_base_controllers)
        end
        
        function [u] = get_control(this,state,time)
            u = rand(1); 
            warning('implement this')
%             control_sequence(time)
        end
        
    end
      
        
end

