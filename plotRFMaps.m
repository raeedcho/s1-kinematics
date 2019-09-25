%% This script plots the receptive field array map, and some other receptive field data

%% set up
    % get the array maps
    if ispc
        dataroot = '';
    else
        dataroot = '/data/raeed';
    end
    mapdir = fullfile(dataroot,'project-data','limblab','s1-kinematics','elec-maps');
    
    rf_map = load(fullfile(mapdir,'rf_map.mat'));
    rf_map = rf_map.rf_map;
    elec_map = load(fullfile(mapdir,'elec-map.mat'));
    elec_map = elec_map.elec_map;
    
    % attach rf color to map
    rf_map = join(rf_map,getRFColorTable());
    
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

%% plots
    % plot the array map
        plotArrayMap(rf_map,struct('map_plot','receptive_field_color'))
    
    % plot out a distalness vs. lateralness plot
        mapping_sessions = unique(array_map(:,{'monkey','date'}));
        [monkeys,~,monkey_idx] = unique(mapping_sessions(:,'monkey'));
        num_sessions = histcounts(monkey_idx);
        session_dates = cell(1,3);
        for monkeynum = 1:height(monkeys)
            session_dates{monkeynum} = mapping_sessions.date(monkey_idx==monkeynum);
        end
        
        % Figure out subplot situation (monkeys on different rows, sessions on different columns)
        num_subplot_cols = max(num_sessions);
        num_subplot_rows = height(monkeys);
        
        % loop over monkeys/days
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:height(monkeys)
            for sessionnum = 1:num_sessions(monkeynum)
                % get map table for given session
                [~,rf_map_session] = getNTidx(rf_map,...
                    'monkey',monkeys{monkeynum,1},...
                    'date',session_dates{monkeynum}{sessionnum});
        
                % set subplot
                plotnum = (monkeynum-1)*num_subplot_cols+sessionnum;
                subplot(num_subplot_rows,num_subplot_cols,plotnum)
        
                % scatter plot the distalness against lateralness
                scatter3(...
                    rf_map_session.elec_anteriorness,...
                    rf_map_session.elec_lateralness,...
                    rf_map_session.rf_distalness,...
                    [],'k','filled')
        
                title(vertcat(monkeys{monkeynum,1},session_dates{monkeynum}(sessionnum)))
                xlabel('Electrode lateralness (au)')
                ylabel('Receptive field distalness (au)')

                set(gca,'box','off','tickdir','out')
            end
        end
