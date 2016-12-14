clear
close all
clc

here = fileparts(mfilename('fullpath'));
beats_folder = fullfile(here,'beats');
addpath(fullfile(here,'BaseBlock'))
addpath(fullfile(here,'beats'))
addpath(fullfile(here,'DataClasses'))
addpath(genpath(fullfile(here,'lib')))
% addpath(fullfile(here,'lib','ppt'))
% addpath(fullfile(here,'lib','xml_io_tools_2007_07'))
% addpath(fullfile(here,'lib','isfieldRecursive'))
% addpath(fullfile(here,'lib','L1_specs','mo','generate_mo'))
addpath(fullfile(here,'ParallelBlock'))
javaaddpath(fullfile(here,'beats','beats-0.1-SNAPSHOT-jar-with-dependencies.jar'))

model_dt = 4;
record_dt = 300;

% create real world model
model = Model(fullfile(beats_folder,'x.xml'),model_dt);

allday_demands = model.get_demands(20000,300,50000);

random_state = model.get_random_state;
zero_state = model.get_zero_state;
num_ctrl_seq=[10 50 40];

rate = 50;
simple_controller = SimpleController(model,rate);
% simple_controller.compute_control_sequence(random_state,allday_demands)

alinea_param.gain_kph = 100;
alinea_param.tgt_dty_vpk = 120;
alinea_param.u_max_vph = 900;
alinea_param.u_min_vph = 100;
alinea = AlineaFeedbackController(model,alinea_param);

[control_sequence,state_trajectory] = model.run_with_controller(alinea,zero_state,allday_demands,record_dt);
 
