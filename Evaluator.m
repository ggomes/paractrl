classdef Evaluator < handle
    
    properties
        controllers
        model
        scores
    end
    
    methods
        
        function this = Evaluator(baseblock,parallelblock,config_file,perturbation,model_dt)
            this.controllers = [baseblock.controllers,parallelblock.controllers];
            this.model = Model(config_file,perturbation,model_dt);
            this.scores = nan(1,length(this.controllers));
        end
        
        function best_controller = choose_best_controller(this,initial_state,predicted_demands)
           
            % evaluate the performance of each controller
            this.scores = nan(1,length(this.controllers));
            for i = 1:length(this.controllers)
                [~,state_trajectory] = this.model.run_with_controller(this.controllers{i},initial_state,predicted_demands);
                this.scores(i) = this.evaluate_performance(state_trajectory);
            end
            
            % pick the best controller
            [~,best_ind] = min(this.scores);
            best_controller = this.controllers{best_ind};
            
        end
                
        function [J] = evaluate_performance(this,state_trajectory)
            
            J = 5;
            warning('implement this')
            
        end

    end
    
end

