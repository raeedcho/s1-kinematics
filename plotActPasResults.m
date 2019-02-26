%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results saved by calculateActPasSeparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up plotting variables
    datadir = '/home/raeed/data/project-data/limblab/s1-kinematics/Results/Separability';
    % filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    files = dir(fullfile(datadir,'*separationResults_run20190224.mat'));
    filename = horzcat({files.name});
    
    % for figure saving
    figdir = '/home/raeed/Wiki/Projects/limblab/s1-kinematics/figures/Encoding';
    run_date = char(datetime('today','format','yyyyMMdd'));

    monkey_names = {'Chips','Han'};
    models_to_plot = {'ext','extforce','handelbow'};

    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Compile information over all files
    [model_eval,model_tuning,tuning_corr,shift_vaf,tuned_neurons] = deal(cell(length(monkey_names),size(session_colors,1)));
    session_ctr = zeros(length(monkey_names),1);
    fileclock = tic;
    fprintf('Started loading files...\n')
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % classify monkey and session number
        monkey_idx = find(strcmpi(encoderResults.crossEval.monkey{1},monkey_names));
        session_ctr(monkey_idx) = session_ctr(monkey_idx) + 1;

        % extract...stuff

        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(filename),toc(fileclock))
    end
