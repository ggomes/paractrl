% main script for generating object factories from xsd files using
% xsd_to_mo

clc
clear
close all

folder_utilities = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(folder_utilities,'lib','xml_io_tools_2007_07'))

folder_model_objects = fileparts(mfilename('fullpath'));
folder_xsd = fullfile(folder_model_objects,'xsd');
folder_output = fullfile(folder_model_objects,'generate_mo');

% schema files .......................................

xsdfile{1}{1} = fullfile(folder_xsd,'scenario_full.xsd');

% % version 2
% xsdfile{2}{1} = fullfile(folder_xsd,'exe-params.xsd');
% xsdfile{2}{2} = fullfile(folder_xsd,'generic.xsd');
% xsdfile{2}{3} = fullfile(folder_xsd,'measurements.xsd');
% xsdfile{2}{4} = fullfile(folder_xsd,'network.xsd');
% xsdfile{2}{5} = fullfile(folder_xsd,'scenario.xsd');
% xsdfile{2}{6} = fullfile(folder_xsd,'traffic-state.xsd');

% run xsd_to_mo on each one ..........................
xsd_to_mo(xsdfile,folder_output,'generate_mo')

disp('done')
