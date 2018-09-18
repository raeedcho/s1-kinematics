function model_colors = getModelColors(model_alias)
% get model colors, given a list of aliases
    if ~iscell(model_alias)
        model_alias = {model_alias};
    end
    model_colors = zeros(length(model_alias),3);
    for modelnum = 1:length(model_alias)
        switch model_alias{modelnum}
        case 'ext'
            model_colors(modelnum,:) = [247, 148, 30]/255;
        case 'opensim_ext'
            model_colors(modelnum,:) = [247, 148, 30]/255;
        case 'ego'
            model_colors(modelnum,:) = [105, 189, 69]/255;
        case 'opensim_ego'
            model_colors(modelnum,:) = [105, 189, 69]/255;
        case 'cyl'
            model_colors(modelnum,:) = [113, 191, 110]/255;
        case 'opensim_cyl'
            model_colors(modelnum,:) = [113, 191, 110]/255;
        case 'joint'
            model_colors(modelnum,:) = [38, 34, 98]/255;
        case 'musc'
            model_colors(modelnum,:) = [0, 174, 239]/255;
        case 'handelbow'
            model_colors(modelnum,:) = [193, 25, 47]/255;
        case 'opensim_handelbow'
            model_colors(modelnum,:) = [179, 44, 224]/255;
        case 'ego_handelbow'
            model_colors(modelnum,:) = [119, 255, 189]/255;
        end
    end

