function plotArrayMaps(map_plot)
% This function plots array modality maps and receptive field maps

%% Set up
% modality colors
modality_colors = [...
    252,141,89;...
    255,255,191;...
    153,213,148]/255;
modalities = {'cutaneous';'mixed';'proprioceptive'};
modality_color_table = table(modalities,modality_colors,'VariableNames',{'modality','modality_color'});

% receptive field colors
rf_colors = [...
    0,0,0;...
    253,224,221;...
    252,197,192;...
    250,159,181;...
    247,104,161;...
    221,52,151;...
    174,1,126;...
    122,1,119]/255;
% rf_colors = [...
%     255,247,243;...
%     253,224,221;...
%     252,197,192;...
%     250,159,181;...
%     247,104,161;...
%     221,52,151;...
%     174,1,126;...
%     122,1,119]/255;
rf_list = {'arm';'torso';'shoulder';'humerus';'elbow';'forearm';'wrist';'hand'};
rf_color_table = table(rf_list,rf_colors,'VariableNames',{'receptive_field','rf_color'});

% array rotations (such that medial is up and anterior is to the right)
array_rot = {...
    'Chips',[-55 90];...
    'Han',[-30 -90];...
    'Lando',[0 -90]};
array_rot = cell2table(array_rot,'VariableNames',{'monkey','array_rotation'});

if ispc
    dataroot = '';
else
    dataroot = '/data/raeed';
end

mapdir = fullfile(dataroot,'project-data','limblab','s1-kinematics','elec-maps');
savedir = fullfile(dataroot,'project-data','limblab','s1-kinematics','Results','arraymaps');

% load in electrode maps
elec_map = load(fullfile(mapdir,'elec-map.mat'));
elec_map = elec_map.elec_map;
rf_map = load(fullfile(mapdir,'rf-map.mat'));
rf_map = rf_map.rf_map;

% attach array locations to the receptive fields
array_map = join(rf_map,elec_map);

% attach modality colors to receptive fields
array_map = join(array_map,modality_color_table);

% attach receptive field colors
array_map = join(array_map,rf_color_table);

% attach array rotations
array_map = join(array_map,array_rot);

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
        if strcmpi(map_plot,'modality')
            color_table = array_map_session.modality_color;
        elseif strcmpi(map_plot,'rf')
            color_table = array_map_session.rf_color;
        else
            error('map_plot must be either ''modality'' or ''rf''')
        end
        for tabrownum = 1:height(array_map_session)
            rectangle(...
                'Position',[array_map_session.colNum(tabrownum) array_map_session.rowNum(tabrownum) 0.8 0.8],...
                'FaceColor',color_table(tabrownum,:))
        end

        title(vertcat(monkeys{monkeynum,1},session_dates{monkeynum}(sessionnum)))
        axis image
        axis off
        view(array_map_session.array_rotation(1,:))
    end
end
subplot(num_subplot_rows,num_subplot_cols,num_subplot_rows*num_subplot_cols)
if strcmpi(map_plot,'modality')
    color_table = modality_color_table;
elseif strcmpi(map_plot,'rf')
    color_table = rf_color_table;
else
    error('map_plot must be either ''modality'' or ''rf''')
end
% color_table = rf_color_table;
for colornum = 1:height(color_table)
    rectangle('Position',[1 colornum 0.8 0.8],'FaceColor',color_table{colornum,2})
    text(2,colornum+0.4,color_table{colornum,1},'FontSize',14)
end
rectangle('Position',[1 0 0.8 0.8],'Curvature',1)
text(2,0.4,'unmapped/none found','FontSize',14)
axis image
axis off
