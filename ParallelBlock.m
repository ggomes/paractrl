classdef ParallelBlock < handle
    %UNTITLED15 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        controllers 
        ctrl_seq_input % array of controlsequence, as large as the number of base block controllers
        warm_start_control_sequence % this is an array of control sequences
                                    % set by set_warm_start
    end
    
    methods
        
        function this = ParallelBlock(controllers)
            this.controllers = controllers;
        end
        
        function this = compute_control_sequence(this,warm_starts,current_state,predicted_demands)
            for i=1:length(this.controllers)
                c = this.controllers{i};
                for j=1:length(warm_starts)
                    c.control_sequence(j) = c.compute_control_sequence(warm_starts(j),current_state,predicted_demands);
                end
            end
        end
            

        
    end
    
end

