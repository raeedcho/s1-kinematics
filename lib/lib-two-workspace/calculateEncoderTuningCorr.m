function tuning_corr = calculateEncoderTuningCorr(encoderResults,params)
    % Calculates correlation between model tuning curves and actual tuning curves
    % Inputs:
    %   encoderResults - The struct output from mwEncoders
    %   params - Paramers struct
    %       .model_aliases - aliases for models names to evaluate
    %       .neural_signal - name of actual neural signal
    % Outputs:
    %   tuning_corr - table with correlation values for each crossval sample
    %       Table will match neuron table standard, with columns for:
    %           monkey (meta)
    %           date (meta)
    %           task (meta)
    %           signalID (meta)
    %           crossvalID (meta) - [repeatnum foldnum]
    %           {model_aliases}_tuningCorr
    %
    % Note: this code assumes that the crossTuning table is in a specific order, i.e.
    % neurons are nested inside spacenum, nested inside folds, nested inside repeats
    % and order inside nestings is preserved and identical

    % set default params
    model_aliases = encoderResults.params.model_aliases;
    neural_signal = 'S1_FR';
    if nargin>1
        assert(isstruct(params),'params should be a struct!')
        assignParams(who,params);
    end

    % define some variables for convenience
    [~,crossTuningPM] = getNTidx(encoderResults.crossTuning,'spaceNum',1);
    [~,crossTuningDL] = getNTidx(encoderResults.crossTuning,'spaceNum',2);

    % figure out what meta columns to keep
    metacols = strcmpi(crossTuningPM.Properties.VariableDescriptions,'meta') &...
        ~strcmpi(crossTuningPM.Properties.VariableNames,'spacenum');

    % loop through tuning table
    tuning_corr_cell = cell(height(crossTuningPM),1);
    for tablerownum = 1:height(crossTuningPM)
        real_tuning_shape = horzcat(...
            crossTuningPM(tablerownum,:).(sprintf('%s_velCurve',neural_signal)),...
            crossTuningDL(tablerownum,:).(sprintf('%s_velCurve',neural_signal)))';

        model_tuning_shape = zeros(length(real_tuning_shape),length(model_aliases));
        for modelnum = 1:length(model_aliases)
            % get tuning shapes
            model_tuning_shape(:,modelnum) = horzcat(...
                crossTuningPM(tablerownum,:).(sprintf('glm_%s_model_velCurve',model_aliases{modelnum})),...
                crossTuningDL(tablerownum,:).(sprintf('glm_%s_model_velCurve',model_aliases{modelnum})))';
        end
        % get correlation value
        covar_mat = nancov([model_tuning_shape real_tuning_shape]);
        model_corr = covar_mat(end,1:end-1)./sqrt(diag(covar_mat(1:end-1,1:end-1))'*covar_mat(end,end));
        model_entry = array2table(model_corr,...
            'VariableNames',strcat(model_aliases,'_tuningCorr'));
        model_entry.Properties.VariableDescriptions = repmat({'linear'},size(model_aliases));

        % assemble table
        tuning_corr_cell{tablerownum} = horzcat(...
            crossTuningPM(tablerownum,metacols),...
            model_entry);
    end

    tuning_corr = vertcat(tuning_corr_cell{:});
