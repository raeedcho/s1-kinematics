function model_titles = getModelTitles(model_alias)
% get model titles, given a list of aliases
    if ~iscell(model_alias)
        model_alias = {model_alias};
    end
    model_titles = cell(length(model_alias),1);
    for modelnum = 1:length(model_alias)
        switch model_alias{modelnum}
        case 'ext'
            model_titles{modelnum} = 'Hand-based';
        case 'opensim_ext'
            model_titles{modelnum} = 'OpenSim Hand-based';
        case 'ego'
            model_titles{modelnum} = 'Egocentric';
        case 'opensim_ego'
            model_titles{modelnum} = 'OpenSim Egocentric';
        case 'cyl'
            model_titles{modelnum} = 'Cylindrical ego';
        case 'opensim_cyl'
            model_titles{modelnum} = 'OpenSim Cylindrical ego';
        case 'joint'
            model_titles{modelnum} = 'Joint-based';
        case 'musc'
            model_titles{modelnum} = 'Muscle-based';
        case 'handelbow'
            model_titles{modelnum} = 'Hand/Elbow-based';
        case 'opensim_handelbow'
            model_titles{modelnum} = 'OpenSim Hand/Elbow-based';
        case 'ego_handelbow'
            model_titles{modelnum} = 'Spherical Hand/Elbow-based';
        end
    end

