function plotArrayMap(array_map,params)
% This function plots array map based on params.map_plot (default: 'modality')

%% Set up
if ispc
    dataroot = '';
else
    dataroot = '/data/raeed';
end
mapdir = fullfile(dataroot,'project-data','limblab','s1-kinematics','elec-maps');

% which map to plot
map_plot = 'modality';
clims = [];

if nargin>1
    assignParams(who,params)
end

% integrety check
assert(any(strcmp(sprintf('%s_color',map_plot),array_map.Properties.VariableNames)))

%% join array map with electrode map
% load in electrode maps
elec_map = load(fullfile(mapdir,'elec-map.mat'));
elec_map = elec_map.elec_map;

% attach array locations to the receptive fields
array_map = join(array_map,elec_map);

% attach array rotations
array_map = join(array_map,getArrayRotationTable());

% get mapping sessions
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
        [~,array_map_session] = getNTidx(array_map,...
            'monkey',monkeys{monkeynum,1},...
            'date',session_dates{monkeynum}{sessionnum});

        % set subplot
        plotnum = (monkeynum-1)*num_subplot_cols+sessionnum;
        subplot(num_subplot_rows,num_subplot_cols,plotnum)

        % first plot all electrode locations for array as boxes
        [~,elec_map_monkey] = getNTidx(elec_map,'monkey',monkeys{monkeynum,1});
        for tabrownum = 1:height(elec_map_monkey)
            rectangle(...
                'Position',[elec_map_monkey.colNum(tabrownum) elec_map_monkey.rowNum(tabrownum) 0.8 0.8],...
                'Curvature',1)
            hold on
        end

        % draw over rectangles with colors based on modality
        color_table = array_map_session.(sprintf('%s_color',map_plot));
        for tabrownum = 1:height(array_map_session)
            rectanglePatch(...
                [array_map_session.colNum(tabrownum) array_map_session.rowNum(tabrownum) 0.8 0.8],...
                color_table(tabrownum,:))
        end

        if ~strcmpi(map_plot,'modality') && ~strcmpi(map_plot,'receptive_field')
            if ~isempty(clims)
                caxis(clims)
            end
            colormap(viridis)
            colorbar
        end

        title(vertcat(monkeys{monkeynum,1},session_dates{monkeynum}(sessionnum)))
        axis image
        axis off
        view(array_map_session.array_rotation(1,:))
    end
end

% print out legend
if strcmpi(map_plot,'modality')
    color_table = getModalityColorTable();
elseif strcmpi(map_plot,'receptive_field')
    color_table = getRFColorTable();
else
    return
    % error('map_plot must be either ''modality'' or ''rf''')
end

subplot(num_subplot_rows,num_subplot_cols,num_subplot_rows*num_subplot_cols)
for colornum = 1:height(color_table)
    rectangle('Position',[1 colornum 0.8 0.8],'FaceColor',color_table{colornum,2})
    text(2,colornum+0.4,color_table{colornum,1},'FontSize',14)
end
rectangle('Position',[1 0 0.8 0.8],'Curvature',1)
text(2,0.4,'unmapped/none found','FontSize',14)
axis image
axis off
