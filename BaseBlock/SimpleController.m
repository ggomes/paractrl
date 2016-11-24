classdef SimpleController < AbstractController
    
    properties
        rate
    end
    
    methods
        
        function this = SimpleController(rate)
            this@AbstractController(1);
            this.rate = rate;
        end
            
        
        function [u] = get_control(this,state,time)
            u = this.rate;
        end
        
    end
    
end

