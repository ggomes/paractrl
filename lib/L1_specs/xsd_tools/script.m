clear 
close all

L0root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
L1 = fullfile(L0root,'L1_specs','mo','xsd','L1_scenario.xsd');
beats = fullfile(L0root,'simulation','beats.xsd');

[result,why] = compare_xsd(L1,beats)