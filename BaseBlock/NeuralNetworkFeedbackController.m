classdef NeuralNetworkFeedbackController < FeedbackController
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function this = NeuralNetworkFeedbackController(model)
            this@FeedbackController(model);
        end
            
        
        function [u] = get_control(this,state,time)
           
            % correct this later
%             u = evaluate_neural_network(state);
            u = rand(1);
            
        end
        
    end
    
end

