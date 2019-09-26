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

    % Test out predicted vs actual distalness
        % figure('defaultaxesfontsize',18)
        % for monkeynum = 1:height(monkeys)
        %     for sessionnum = 1:num_sessions(monkeynum)
        %         % get current session
        %         [~,rf_map_session] = getNTidx(rf_map,...
        %             'monkey',monkeys{monkeynum,1},...
        %             'date',session_dates{monkeynum}{sessionnum});

        %         % set subplot
        %         plotnum = (monkeynum-1)*num_subplot_cols+sessionnum;
        %         subplot(num_subplot_rows,num_subplot_cols,plotnum)

        %         lm = rf_map_session.distalness_model{1};
        %         scatter(...
        %             lm.predict(rf_map_session(:,{'rowNum','colNum'})),...
        %             rf_map_session.rf_distalness,...
        %             [],'k','filled')
        %         hold on
        %         plot([0 7],[0 7],'--k','linewidth',2)
        %         title(vertcat(monkeys{monkeynum,1},session_dates{monkeynum}(sessionnum)))
        %         xlabel('Predicted rf distalness (au)')
        %         ylabel('Receptive field distalness (au)')
        %         set(gca,'box','off','tickdir','out',...
        %             'xlim',[0 7],'ylim',[0 7])
        %     end
        % end

