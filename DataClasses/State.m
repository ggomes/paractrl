classdef State
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        link_ids    % I x 1
        link_veh    % I x 1 [veh / link]
    end
    
    methods( Access=public ) 
        
        function this = State(link_ids)
            this.link_ids = link_ids;
            if size(this.link_ids,1)==1
                this.link_ids = this.link_ids';
            end
            this.link_veh = [];
        end
        
        function [x]=get_veh_for_links(this,link_ids)
            ind = index_into(link_ids,this.link_ids);
            if isempty(ind)
                error('link id not present')
            end
            x = this.link_veh(ind);
        end
        
    end
    
end

