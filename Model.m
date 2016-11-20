classdef Model < handle
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dt                      % seconds
        current_time;           % seconds after midnight
        link_ids                % 1 x I ... link ids
        controlled_link_ids     % 1 x Ic ... list of links that have ramp meters on them
        source_link_ids         % 1 x Is ... list o links that are sources of demand. 
        demands                 % structure array with demand information ordered like source_link_ids
        beats
    end
    
    methods(Access=public)

        function this = Model(config_file,dt)
            this.dt = dt;

            % load the beats model
            this.beats = BeatsSimulation;
            this.beats.load_scenario(config_file);
            this.beats.scenario_ptr.change_units_to('SI');
            this.beats.create_beats_object(struct('SIM_DT',dt))
            this.beats.intialize_beats_persistent(struct('SIM_DT',4,'OUTPUT_DT',300));

            this.link_ids = this.beats.scenario_ptr.get_link_ids;
            
            if ~isfield(this.beats.scenario_ptr.scenario,'ActuatorSet')
                error('1')
            end
            
            if ~isfield(this.beats.scenario_ptr.scenario.ActuatorSet,'actuator')
                error('2')
            end
            
            % read ramp meters
            actuator_link_ids = nan(2,length(this.beats.scenario_ptr.scenario.ActuatorSet.actuator));
            for i=1:length(this.beats.scenario_ptr.scenario.ActuatorSet.actuator)
                act = this.beats.scenario_ptr.scenario.ActuatorSet.actuator(i);
                if ~strcmp(act.actuator_type.ATTRIBUTE.name,'ramp_meter')
                    continue
                end
                actuator_link_ids(1,i) = act.ATTRIBUTE.id;
                actuator_link_ids(2,i) = act.scenarioElement.ATTRIBUTE.id;
            end
            
            if ~isfield(this.beats.scenario_ptr.scenario,'ControllerSet')
                error('3')
            end
            
            if ~isfield(this.beats.scenario_ptr.scenario.ControllerSet,'controller')
                error('4')
            end
            
            if numel(this.beats.scenario_ptr.scenario.ControllerSet.controller)~=1
                error('5')
            end
            
            if numel(this.beats.scenario_ptr.scenario.ControllerSet.controller.target_actuators)~=1
                error('5')
            end
            
            ctrl_act_ids = arrayfun(@(x) x.ATTRIBUTE.id ,this.beats.scenario_ptr.scenario.ControllerSet.controller.target_actuators.target_actuator);
            this.controlled_link_ids = actuator_link_ids(2,ismember(actuator_link_ids(1,:),ctrl_act_ids));
            this.source_link_ids =this.link_ids(this.beats.scenario_ptr.is_source_link);
            
            % store demand profiles
            demand_prof = this.beats.scenario_ptr.get_demandprofiles_with_IDs(this.source_link_ids);
            numdemand = numel(demand_prof);
            
            this.demands = repmat(struct('link_id',nan,'time',[],'vph',[]),1,numel(this.source_link_ids));
            for i=1:numdemand
                ind = this.source_link_ids==demand_prof(i).ATTRIBUTE.link_id_org;
                this.demands(ind).link_id = demand_prof(i).ATTRIBUTE.link_id_org;
                this.demands(ind).vph = demand_prof(i).demand.CONTENT;
                endtime = demand_prof(i).ATTRIBUTE.start_time + demand_prof(i).ATTRIBUTE.dt*(numel(this.demands(ind).vph)-1);
                this.demands(ind).time = demand_prof(i).ATTRIBUTE.start_time : demand_prof(i).ATTRIBUTE.dt : endtime;
            end

        end
        
        % pass in:
        %    controller @AbstractController 
        %    initial_state @State 
        %    predicted_demands @DemandProfile 
        % returns :
        %    control_sequence @ControlSequence
        %    state_trajectory @StateTrajectory 
        function [control_sequence,state_trajectory]=run_with_controller(this,controller,initial_state,predicted_demands,record_dt)
            
            initial_time = predicted_demands.get_initial_time;
            final_time = predicted_demands.get_final_time;
            
            % initialize the model
            this.set_state(initial_state);
            this.set_time(initial_time);
            
            % run the model            
            state_trajectory = StateTrajectory(this.link_ids);
            
            if isobject(controller)
                control_sequence = ControlSequence(controller,this.controlled_link_ids);
            else
                control_sequence = nan;
            end
            
            while this.current_time < final_time
                    
                if mod(this.current_time,record_dt)==0
                    fprintf('%f\t%f\n',this.current_time,final_time)
                end
                
                % get and record the state
                current_state = this.get_state;
                
                if mod(this.current_time,record_dt)==0
                    state_trajectory.add_values(this.current_time,current_state);
                end
                
                % get and record the control
                if isobject(controller)
                    current_control = controller.get_control(current_state,this.current_time);
                    this.set_control(current_control);
                    if mod(this.current_time,record_dt)==0
                        control_sequence.add_values(this.current_time,current_control);
                    end
                end

                % advance the simulation
                this.advance;
            end
                        
        end
        
        % pass in:
        %    start_time         [seconds after midnight]
        %    end_time           [seconds after midnight]
        % returns :
        %    x @DemandProfile
        function [x] = get_demands(this,start_time,sample_dt,end_time)   
            x = DemandProfile(this.source_link_ids);
            x.time = start_time:sample_dt:end_time;
            x.demands = nan(numel(this.source_link_ids),numel(x.time));
            for i=1:length(this.source_link_ids)
                x.demands(i,:) = interp1(this.demands(i).time,this.demands(i).vph,x.time,'previous','extrap');
            end
        end
        
        % returns a State object with random (but feasible) values
        % can be used for debugging
        function x = get_random_state(this)
            x = State(this.link_ids);
            fds = this.beats.scenario_ptr.get_fds_with_linkIDs(this.link_ids);
            link_lanes = this.beats.scenario_ptr.get_link_lanes;
            link_lengths = this.beats.scenario_ptr.get_link_lengths;
            max_veh = link_lanes.*link_lengths.*[fds.capacity].*( 1./[fds.congestion_speed] + 1./[fds.free_flow_speed] );
            x.link_veh = rand(1,numel(this.link_ids)).*max_veh;
        end
        
        % returns a State object with no vehicles
        function x = get_zero_state(this)
            x = State(this.link_ids);
            x.link_veh = zeros(1,numel(this.link_ids));
        end
        
    end
    
    methods(Access=private)
                
        function []=set_state(this,initial_state)            
            this.beats.beats.set.totalDensity(initial_state.link_veh');
        end
        
        function []=set_time(this,t)
            this.current_time = t;
        end
        
        function [x]=get_state(this)
            x = State(this.link_ids);
            x.link_veh = this.beats.beats.get.totalDensity(nan);
        end
        
        function []=set_control(this,u)
            warning('implement this')
        end
        
        function []=advance(this)
            this.beats.beats.advanceNSeconds(this.dt);
            this.current_time = this.current_time + this.dt;
        end
 
    end
    
end

