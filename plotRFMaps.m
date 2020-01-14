%% This script plots the receptive field array map, and some other receptive field data

%% set up
    % get the array maps
    if ispc
        dataroot = 'G:\raeed\project-data\limblab\s1-kinematics';
    else
        dataroot = '/data/raeed/project-data/limblab/s1-kinematics';
    end
    mapdir = fullfile(dataroot,'sensory_mapping');

    rf_map = readtable(fullfile(mapdir,'rf_map.csv'));
    elec_map = readtable(fullfile(mapdir,'elec_map.csv'));
    
    % attach rf color to map
    rf_map = join(rf_map,getRFColorTable());

    % attach modality color to map
    rf_map = join(rf_map,getModalityColorTable());
    
    % attach rf distalness to map
    rf_map = join(rf_map,getRFDistalnessTable());

    % attach electrode table
    rf_map = join(rf_map,elec_map);
    
%% array plots
    % plot the array modality map
        figure('defaultaxesfontsize',18)
        plotArrayMap(rf_map,struct('mapdir',mapdir,'map_plot','modality_color'));

    % plot out the rf map
        rf_colormap = getRFColorTable();
        figure('defaultaxesfontsize',18)
        lm_table = plotArrayMap(rf_map,struct(...
            'mapdir',mapdir,...
            'map_plot','rf_distalness',...
            'clims',[0 7],...
            'calc_linmodels',true,...
            'show_colorbar',false,...
            'cmap',rf_colormap.receptive_field_color));

