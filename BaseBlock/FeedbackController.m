classdef FeedbackController < AbstractController
    
    properties
        model   % model used to generate a control sequence
    end
    
    methods
        
        function this = FeedbackController(model)
            this@AbstractController(1);
            this.model = model;
        end
        
        function [] = compute_control_sequence(this,initial_state,predicted_demands)
            this.control_sequence = this.model.run_with_controller(this,initial_state,predicted_demands);   
        end
        
    end
    
end

