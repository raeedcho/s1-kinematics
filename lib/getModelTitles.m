function model_titles = getModelColors(model_alias)
% get model titles, given a list of aliases
    if ~iscell(model_alias)
        model_alias = {model_alias};
    end
    model_titles = cell(length(model_alias),1);
    for modelnum = 1:length(model_alias)
        switch model_alias{modelnum}
        case 'musc'
            model_titles{modelnum} = 'Muscle-based';
        case 'ext'
            model_titles{modelnum} = 'Hand-based';
        case 'ego'
            model_titles{modelnum} = 'Egocentric';
        case 'cyl'
            model_titles{modelnum} = 'Cylindrical ego';
        case 'joint'
            model_titles{modelnum} = 'Joint-based';
        case 'markers'
            model_titles{modelnum} = 'Hand/Elbow-based';
        end
    end

