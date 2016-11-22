clear
close all

here = fileparts(mfilename('fullpath'));
beats_folder = fullfile(here,'beats');
addpath(fullfile(here,'BaseBlock'))
addpath(fullfile(here,'beats'))
addpath(fullfile(here,'DataClasses'))
addpath(fullfile(here,'lib','ppt'))
addpath(fullfile(here,'lib','xml_io_tools_2007_07'))
addpath(fullfile(here,'ParallelBlock'))
javaaddpath(fullfile(here,'beats','beats-0.1-SNAPSHOT-jar-with-dependencies.jar'))

model_dt = 4;
record_dt = 300;

% create real world model
model = Model(fullfile(beats_folder,'x.xml'),model_dt);

allday_demands = model.get_demands(20000,300,50000);

random_state = model.get_random_state;
zero_state = model.get_zero_state;

[control_sequence,state_trajectory] = model.run_with_controller(nan,zero_state,allday_demands,record_dt);
 
