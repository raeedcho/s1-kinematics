%% This script plots the receptive field array map, and some other receptive field data

%% set up
    % get the array maps
    if ispc
        dataroot = '';
    else
        dataroot = '/data/raeed';
    end
    mapdir = fullfile(dataroot,'project-data','limblab','s1-kinematics','elec-maps');
    
    rf_map = load(fullfile(mapdir,'rf-map.mat'));
    rf_map = rf_map.rf_map;
    elec_map = load(fullfile(mapdir,'elec-map.mat'));
    elec_map = elec_map.elec_map;
    
    % attach rf color to map
    rf_map = join(rf_map,getRFColorTable());

    % attach modality color to map
    rf_map = join(rf_map,getModalityColorTable());
    
    % attach rf distalness to map
    rf_map = join(rf_map,getRFDistalnessTable());
    
    % attach array rotation
    rf_map = join(rf_map,getArrayRotationTable());
    
    % attach electrode map
    rf_map = join(rf_map,elec_map);
    
    % attach lateralness information to each channel
    lateralness = -rf_map.colNum.*sind(rf_map.array_rotation(:,1)) + rf_map.rowNum.*cosd(rf_map.array_rotation(:,1));
    rf_map = horzcat(rf_map,table(lateralness,'VariableNames',{'elec_lateralness'}));
    anteriorness = rf_map.colNum.*cosd(rf_map.array_rotation(:,1)) + rf_map.rowNum.*sind(rf_map.array_rotation(:,1));
    rf_map = horzcat(rf_map,table(anteriorness,'VariableNames',{'elec_anteriorness'}));

%% array plots
    % plot the array modality map
        figure('defaultaxesfontsize',18)
        plotArrayMap(rf_map,struct('map_plot','modality_color'));

    % plot out the rf map
        rf_colormap = getRFColorTable();
        figure('defaultaxesfontsize',18)
        lm_table = plotArrayMap(rf_map,struct(...
            'map_plot','rf_distalness',...
            'clims',[0 7],...
            'calc_linmodels',true,...
            'show_colorbar',false,...
            'cmap',rf_colormap.receptive_field_color));

        rf_map = join(rf_map,lm_table);
