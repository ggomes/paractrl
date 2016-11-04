classdef AlineaFeedbackController < FeedbackController
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        gain
    end
    
    methods
        
        function this = AlineaFeedbackController(model,gain)
            this@FeedbackController(model);
            this.gain = gain;
        end
            
        
        function [u] = get_control(this,state,time)
           
            % correct this later
            % u = this.gain*state;
            u = rand(1);
            
        end
        
    end
    
end

