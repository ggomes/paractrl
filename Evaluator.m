classdef Evaluator < handle
    
    properties
        controllers
        model
        best_controller
    end
    
    methods
        
        function this = Evaluator(baseblock,parallelblock,config_file,perturbation,model_dt)
            this.controllers = [baseblock.controllers,parallelblock.controllers];
            this.model = Model(config_file,perturbation,model_dt);
            this.best_controller = nan;
        end
        
        function this = evaluate(this)
           
            
            for i = 1:length(this.controllers)
                this.controllers(i)
            end
            
            
        end
        
        
        
    end
    
end

