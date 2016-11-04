classdef Model < handle
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dt
        current_time;
        beats
    end
    
    methods

        function this = Model(config_file,perturbation,dt)
            this.dt = dt;
            warning('implement this')
        end
        
        
        function [control_sequence,state_trajectory]=run_with_controller(this,controller,initial_state,predicted_demands)
            
            initial_time = predicted_demands.get_initial_time;
            final_time = predicted_demands.get_final_time;
            
            % initialize the model
            this.set_state(initial_state);
            this.set_time(initial_time);
            
            % run the model
            u_profile = [];
            x_profile = initial_state;
            while this.get_time < final_time
                current_state = this.get_state;
                x_profile = [x_profile current_state];
                u_current = controller.get_control(current_state,this.get_time);
                u_profile = [u_profile u_current];
                this.set_control(u_current);
                this.advance;
            end
            
            control_sequence = ControlSequence(controller,u_profile);
            
        end
        
        function []=set_state(this,initial_state)
            warning('implement this')
        end
        
        function []=set_time(this,t)
            this.current_time = t;
        end
        
        function [x]=get_time(this)
            x = this.current_time;
        end
        
        function [x]=get_state(this)
            x = nan;
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
            time = start_time:this.dt:end_time;
            x = DemandProfile(time);
            x.add_values(1,rand(1,length(time)));
            x.add_values(2,rand(1,length(time)));
            warning('implement this')
        end
        
    end
    
end

