function model_titles = getModelTitles(model_alias)
% get model titles, given a list of aliases
    if iscell(model_alias)
        model_alias_cell = model_alias;
    else
        model_alias_cell = {model_alias};
    end
    model_titles = cell(length(model_alias_cell),1);
    for modelnum = 1:length(model_alias_cell)
        switch model_alias_cell{modelnum}
        case 'ext'
            model_titles{modelnum} = 'Hand-only';
        case 'extforce'
            model_titles{modelnum} = 'Hand Kinematic-Force';
        case 'opensim_ext'
            model_titles{modelnum} = 'OpenSim Extrinsic';
        case 'ego'
            model_titles{modelnum} = 'Egocentric';
        case 'opensim_ego'
            model_titles{modelnum} = 'OpenSim Egocentric';
        case 'cyl'
            model_titles{modelnum} = 'Cylindrical ego';
        case 'opensim_cyl'
            model_titles{modelnum} = 'OpenSim Cylindrical ego';
        case 'joint'
            model_titles{modelnum} = 'Joint';
        case 'musc'
            model_titles{modelnum} = 'Muscle';
        case 'handelbow'
            model_titles{modelnum} = 'Whole-arm';
        case 'elbow'
            model_titles{modelnum} = 'Elbow';
        case 'opensim_handelbow'
            model_titles{modelnum} = 'OpenSim Hand-Elbow';
        case 'ego_handelbow'
            model_titles{modelnum} = 'Spherical Hand-Elbow';
        case 'handelbow_surprise'
            model_titles{modelnum} = 'Surprise!';
        end
    end

    if ~iscell(model_alias)
        model_titles = model_titles{1};
    end
