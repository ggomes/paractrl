clear
close all

config_file = '';
evaluator_pert = 0.2;
start_time = 0;
controller_dt = 300;
model_dt = 5;
time_horizon = 300;

% 0. create real world model
model = Model(config_file,0,model_dt);

% 1. Create base block
baseblock = BaseBlock({ ...
    AlineaFeedbackController(model,5) , ...
    NeuralNetworkFeedbackController(model) , ...
});

% 2. create the base block and Parallel block
size_base_block = baseblock.num_controllers;
parallelblock = ParallelBlock({ ...
    LinearizedMPCController(size_base_block) , ...
    ParametrizedMPCController(size_base_block) ...
});
    
% create the evaluator
evaluator = Evaluator(baseblock,parallelblock,config_file,evaluator_pert,model_dt);

% this is within a time loop....
current_time = start_time;

% get current state and predicted demands
current_state = model.get_state;
predicted_demands = model.predict_demands(current_time,current_time+time_horizon);
predicted_demands.perturb;

warm_starts = baseblock.compute_control_sequence(current_state,predicted_demands);

parallelblock.compute_control_sequence(warm_starts,current_state,predicted_demands);


evaluatorblock.evaluate

best_controller = evaluatorblock.get_best_controller

model.set_controller(best_controller)



% BaseBlock = []; % set of feedback and predictive controllers
% ParallelBlock = []; % set of predictivewarmstart controllers

