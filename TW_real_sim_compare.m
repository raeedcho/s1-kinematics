%% Script to compare PD changes from real neurons to Monte Carlo simulated neurons given intrinsic coordinates

%% Han 20160322
folder = '/home/raeed/Projects/limblab/FSMRes/limblab/User_folders/Raeed/Arm Model/Data/Han/experiment_20160322_TW_RW/';
clear options
options.prefix = 'Han_20160322_RW';
% options.opensim_prefix = '';
options.labnum = 6;
% options.dual_array = 0;
% options.time_window = [0 320];
function_name = 'iris_predict';

dbstop if error
output_data_pred = run_data_processing(function_name,folder,options);
dbclear if error