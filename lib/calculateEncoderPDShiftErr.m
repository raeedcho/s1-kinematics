function err = calculateEncoderPDShiftErr(encoderResults)
% Calculates error in predicting preferred direction shifts for given models
    model_aliases = encoderResults.params.model_aliases;
    err_arr = zeros(100,length(model_aliases));
    shift_tables = calculatePDShiftTables(encoderResults);
    for modelnum = 1:length(model_aliases)
        [~,real_shifts] = getNTidx(shift_tables{end},'signalID',encoderResults.tunedNeurons);
        [~,model_shifts] = getNTidx(shift_tables{modelnum},'signalID',encoderResults.tunedNeurons);
        err_arr_all = model_shifts.velPD-real_shifts.velPD;
        for i = 1:100
            err_idx = 1:length(encoderResults.tunedNeurons);
            err_idx = err_idx + (i-1)*length(encoderResults.tunedNeurons);
            % use a 1-cos style error because of circular data
            % This value will range between 0 and 2
            err_arr(i,modelnum) = mean(1-cos(err_arr_all(err_idx)));
        end
    end

    err = array2table(err_arr,'VariableNames',model_aliases);
        
