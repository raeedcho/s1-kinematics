function lm_table = plotArrayMap(array_map,params)
% This function plots array map based on params.map_plot (default: 'modality')

    %% Set up
    if ispc
        dataroot = 'G:\raeed\project-data\limblab\s1-kinematics';
    else
        dataroot = '/data/raeed/project-data/limblab/s1-kinematics';
    end
    
    % parameters we can change with params
    mapdir = fullfile(dataroot,'sensory_mapping'); % directory of electrode map
    map_plot = 'modality_color'; % what to map
    coarseness = 1; % coarseness of map (averages over NxN area of array)
    clims = []; % color limits
    calc_linmodels = false; % whether to calculate and show linear models of map
    cmap = viridis; % color map to use
    show_colorbar = true; % whether to show colorbar legend
    
    if nargin>1
        assignParams(who,params)
    end
    
    % integrety check
    % assert(any(strcmp(sprintf('%s_color',map_plot),array_map.Properties.VariableNames)))
    assert(isnumeric(array_map.(map_plot)),'Column to plot must be numeric (either index into colormap or n x 3 color array)')
    assert(mod(coarseness,1)==0 && coarseness>0, 'coarseness must be an integer greater than 0')
    
    %% join array map with electrode map
    % load in electrode maps
    elec_map = readtable(fullfile(mapdir,'elec_map.csv'));
    
    % attach array locations to the receptive fields
    array_map = join(array_map,elec_map);
    
    % coarsify array map if necessary
    if coarseness>1
        array_map.rowNum = floor((array_map.rowNum-1)/coarseness)+1;
        array_map.colNum = floor((array_map.colNum-1)/coarseness)+1;
        keycols = ismember(array_map.Properties.VariableNames,{'monkey','date','rowNum','colNum'});
        map_col = strcmpi(array_map.Properties.VariableNames,map_plot);
        array_map = neuronAverage(array_map,struct(...
            'keycols',keycols,...
            'ignorecols',~keycols & ~map_col,...
            'do_ci',false));
    end
    
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
    lm_table = cell(max(num_sessions),height(monkeys));
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
            if coarseness==1
                for tabrownum = 1:height(elec_map_monkey)
                    rectangle(...
                        'Position',[elec_map_monkey.colNum(tabrownum) elec_map_monkey.rowNum(tabrownum) 0.8 0.8],...
                        'Curvature',1)
                    hold on
                end
            else
                for colnum = 1:((max(elec_map_monkey.colNum)-1)/coarseness+1)
                    for rownum = 1:((max(elec_map_monkey.rowNum)-1)/coarseness+1)
                        rectangle(...
                            'Position',[colnum rownum 0.8 0.8],...
                            'Curvature',1)
                        hold on
                    end
                end
            end
    
            % draw over rectangles with colors based on modality
            color_table = array_map_session.(map_plot);
            for tabrownum = 1:height(array_map_session)
                rectanglePatch(...
                    [array_map_session.colNum(tabrownum) array_map_session.rowNum(tabrownum) 0.8 0.8],...
                    color_table(tabrownum,:))
            end
    
            if ~startsWith(map_plot,'modality') && ~startsWith(map_plot,'receptive_field')
                if ~isempty(clims)
                    caxis(clims)
                end
                colormap(cmap)
                if show_colorbar
                    colorbar
                end
            end
    
            if calc_linmodels
                % get linear model of distalness against array position and put in a table
                lm = fitlm(array_map_session(:,{'colNum','rowNum',map_plot}));
    
                % plot if lm is good
                if coefTest(lm)<0.05
                    midpoint = [...
                        (max(elec_map_monkey.colNum)-1)/coarseness+1,...
                        (max(elec_map_monkey.rowNum)-1)/coarseness+1]/2+0.9
                    grad = lm.Coefficients.Estimate(2:end);
                    grad = (5/coarseness)*grad/norm(grad);
                    plot(...
                        [midpoint(1) midpoint(1)+grad(1)],...
                        [midpoint(2) midpoint(2)+grad(2)],...
                        'k','linewidth',2)
                end
    
                % save table
                lm_table{sessionnum,monkeynum} = table(...
                    monkeys{monkeynum,1},...
                    session_dates{monkeynum}(sessionnum),...
                    {lm},...
                    'VariableNames',{'monkey','date',strcat(map_plot,'_linmodel')});
            end
    
            title(vertcat(monkeys{monkeynum,1},session_dates{monkeynum}(sessionnum)))
            axis image
            axis off
            view(array_map_session.array_rotation(1,:))
        end
    end
    lm_table = vertcat(lm_table{:});
    
    % print out legend
    if startsWith(map_plot,'modality')
        color_table = getModalityColorTable();
    elseif startsWith(map_plot,'receptive_field')
        color_table = getRFColorTable();
    elseif startsWith(map_plot,'rf_distalness')
        color_table = getRFDistalnessTable();
    else
        return
        % error('map_plot must be either ''modality'' or ''rf''')
    end
    
    subplot(num_subplot_rows,num_subplot_cols,num_subplot_rows*num_subplot_cols)
    for colornum = 1:height(color_table)
        rectanglePatch([1 colornum 0.8 0.8],color_table{colornum,2})
        text(2,colornum+0.4,color_table{colornum,1},'FontSize',14)
    end
    rectangle('Position',[1 0 0.8 0.8],'Curvature',1)
    text(2,0.4,'unmapped/none found','FontSize',14)
    
    if ~isempty(clims)
        caxis(clims)
    end
    axis image
    axis off
end

function rectanglePatch(pos,facecolor)
% This function plots a rectangle with face color using patch (so we can use the colormap if we want)
    x_verts = [pos(1) pos(1)+pos(3) pos(1)+pos(3) pos(1)];
    y_verts = [pos(2) pos(2) pos(2)+pos(4) pos(2)+pos(4)];

    patch(x_verts,y_verts,facecolor)
end

function array_rot = getArrayRotationTable()
% Function to get a table of rotations for each monkey's Left S1 array
    
    % array rotations (such that medial is down and anterior is to the right, i.e. for Left S1)
    array_rot = {...
        'C',[-55 90];...
        'H',[-30 90];...
        'L',[0 90]};
    array_rot = cell2table(array_rot,'VariableNames',{'monkey','array_rotation'});
end
