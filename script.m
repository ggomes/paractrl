clear
close all

here = fileparts(mfilename('fullpath'));
beats_folder = fullfile(here,'beats');

start_time = 0;
controller_dt = 300;
model_dt = 4;
time_horizon = 300;

% create real world model
model = Model(fullfile(beats_folder,'x.xml'),model_dt);

% Create base block
baseblock = BaseBlock({ ...
    AlineaFeedbackController(model,5) , ...
    NeuralNetworkFeedbackController(model) , ...
});

% create the base block and Parallel block
size_base_block = baseblock.num_controllers;
parallelblock = ParallelBlock({ ...
    LinearizedMPCController(size_base_block) , ...
    ParametrizedMPCController(size_base_block) ...
});
    
% create the evaluator
evaluator = Evaluator(baseblock,parallelblock,config_file,model_dt);

% this is within a time loop....
current_time = start_time;

% get current state and predicted demands
current_state = model.get_state;
predicted_demands = model.get_demands(current_time,current_time+time_horizon);
predicted_demands.perturb;

warm_starts = baseblock.compute_control_sequence(current_state,predicted_demands);

parallelblock.compute_control_sequence(warm_starts,current_state,predicted_demands);


best_controller = evaluator.choose_best_controller(current_state,predicted_demands);

% model.set_controller(best_controller)

