function shift_tables = calculatePDShiftTables(encoderResults,model_names)
% calculates PD shifts for given encoder results and given model names
    if nargin<2
        model_names = encoderResults.params.model_names;
    end
    if ~iscell(model_names)
        model_names = {model_names};
    end
    shift_tables = cell(1,length(model_names));
    for modelnum = 1:length(model_names)
        % select tables for each space
        [~,pm_tuningTable] = getNTidx(encoderResults.crossTuning,'spaceNum',1);
        [~,dl_tuningTable] = getNTidx(encoderResults.crossTuning,'spaceNum',2);
        
        % compose shift table for this model/bootstrap sample
        key_cols = contains(pm_tuningTable.Properties.VariableDescriptions,'meta');
        shift_tables{modelnum} = pm_tuningTable(:,key_cols);
        
        % remove spaceNum from columns
        shift_tables{modelnum}.spaceNum = [];
        
        % get PDs from pm and dl
        pm_PDs = pm_tuningTable.([model_names{modelnum} '_velPD']);
        dl_PDs = dl_tuningTable.([model_names{modelnum} '_velPD']);
        dPDs = minusPi2Pi(dl_PDs-pm_PDs);
        
        tab_append = table(dPDs,'VariableNames',{'velPD'});
        tab_append.Properties.VariableDescriptions = {'circular'};
        shift_tables{modelnum} = [shift_tables{modelnum} tab_append];
        
    end
