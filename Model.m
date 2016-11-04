classdef Model < handle
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dt
        current_time;
        link_ids
        controlled_link_ids
        source_link_ids
        beats
    end
    
    methods

        function this = Model(config_file,perturbation,dt)
            this.dt = dt;
            this.link_ids = [1 2 3]';
            this.controlled_link_ids = [1 2]';
            this.source_link_ids = [1 2]';
            warning('implement this')
        end
        
        
        function [control_sequence,state_trajectory]=run_with_controller(this,controller,initial_state,predicted_demands)
            
            initial_time = predicted_demands.get_initial_time;
            final_time = predicted_demands.get_final_time;
            
            % initialize the model
            this.set_state(initial_state);
            this.set_time(initial_time);
            
            % run the model            
            state_trajectory = StateTrajectory(this.link_ids);
            control_sequence = ControlSequence(controller,this.controlled_link_ids);
            
            while this.current_time < final_time
                current_state = this.get_state;
                current_control = controller.get_control(current_state,this.current_time);
                
                state_trajectory.add_values(this.current_time,current_state);
                control_sequence.add_values(this.current_time,current_control);

                this.set_control(current_control);
                this.advance;
            end
                        
        end
        
        function []=set_state(this,initial_state)
            warning('implement this')
        end
        
        function []=set_time(this,t)
            this.current_time = t;
        end
        
        function [x]=get_state(this)
            x = rand(length(this.link_ids),1);
            warning('implement this')
        end
        
        function []=set_control(this,u)
            warning('implement this')
        end
        
        function []=advance(this)
            this.current_time = this.current_time + this.dt;
            warning('implement this')
        end
        
        % returns a DemandProfile
        function [x] = predict_demands(this,start_time,end_time)
            
            x = DemandProfile(this.source_link_ids);
            time = start_time:this.dt:end_time;
            for t = 1:length(time)
                v = rand(length(this.source_link_ids),1);
                x.add_values(time(t),v);
            end
            warning('implement this')
        end
        
    end
    
end

