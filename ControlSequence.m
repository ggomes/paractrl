classdef ControlSequence < handle
    %UNTITLED16 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        controller
        score
        control_sequence
    end
    
    methods
        
        function this = ControlSequence(controller,control_sequence)
            this.controller = controller;
            this.score = nan;
            this.control_sequence = control_sequence;
        end
        
    end
    
end

