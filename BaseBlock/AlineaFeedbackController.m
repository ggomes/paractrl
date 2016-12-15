classdef AlineaFeedbackController < FeedbackController
    
    properties
        u_prev_vph          % Ix1 [veh/hr]
        ml_link_ids         % Ix1 ids for mainline links corresponding to each onramp
        ml_link_length_km   % Ix1 [km] length of the mainline link
        gain_kph            % [Ix1] [km/hr]
        tgt_dty_vpk         % [Ix1] [veh/km]
        u_max_vph           % [Ix1] [vph]
        u_min_vph           % [Ix1] [vph]
    end
    
    methods( Access=public )
        
        function this = AlineaFeedbackController(model,p)
            this@FeedbackController(model);
            
            % number of controlled onramps
            n = numel(this.controlled_link_ids);
            
            % create parameters for each controlled ramp
            % for now just assign parameters uniformly
            % we may consider in the future having different parameters 
            % for different ramps. 
            this.gain_kph = p.gain_kph*ones(n,1);
            this.tgt_dty_vpk = p.tgt_dty_vpk*ones(n,1);
            this.u_max_vph = p.u_max_vph*ones(n,1);
            this.u_min_vph = p.u_min_vph*ones(n,1);
            
            % find the mainline links corresponding to each onramp
            % this part makes use of BeATS functionality
            link_id_begin_end = model.beats.scenario_ptr.get_link_id_begin_end;
            this.ml_link_ids = nan(n,1);
            this.ml_link_length_km = nan(n,1);
            all_link_length_km = model.beats.scenario_ptr.get_link_lengths('si')/1000;
            for i=1:n
                or_link = this.controlled_link_ids(i);
                merge_node = link_id_begin_end(link_id_begin_end(:,1)==or_link,3);
                up_links_ind = link_id_begin_end(:,3)==merge_node;
                or_ind = model.link_ids==or_link;
                up_ml_link_ind = up_links_ind & ~or_ind;
                if sum(up_ml_link_ind)~=1
                    error('ramp does not have a unique upstream mainline link')
                end
                this.ml_link_ids(i) = model.link_ids(up_ml_link_ind);
                this.ml_link_length_km(i) = all_link_length_km(up_ml_link_ind);
            end
            
            % initialize u to the maximum rate
            this.u_prev_vph = nan(n,1);
            for i=1:n
                this.u_prev_vph(i) = this.u_max_vph(i);
            end
            
        end
            
        % returns a singleton ControlSequence
        function [u] = get_control(this,state,time)
                        
            % extract the mainline densities from the state
            ml_veh = state.get_veh_for_links(this.ml_link_ids);
            ml_dty_vpk = ml_veh./this.ml_link_length_km;
                
            % alinea control law
            u_vph = this.u_prev_vph + this.gain_kph.*(this.tgt_dty_vpk - ml_dty_vpk);
            
            % max and min
            u_vph = min([u_vph this.u_max_vph],[],2);
            u_vph = max([u_vph this.u_min_vph],[],2);
                
            % record u_prev_vph
            this.u_prev_vph = u_vph;
                
            % cast to a ControlSequence
            u = ControlSequence(this,this.controlled_link_ids);
            u.add_values(time,u_vph);
        end
        
    end
    
end

