function model_colors = getModelColors(model_alias)
% get model colors, given a list of aliases
    if ~iscell(model_alias)
        model_alias = {model_alias};
    end
    model_colors = zeros(length(model_alias),3);
    for modelnum = 1:length(model_alias)
        switch model_alias{modelnum}
        case 'musc'
            model_colors(modelnum,:) = [0, 174, 239]/255;
        case 'ext'
            model_colors(modelnum,:) = [247, 148, 30]/255;
        case 'ego'
            model_colors(modelnum,:) = [105, 189, 69]/255;
        case 'cyl'
            model_colors(modelnum,:) = [113, 191, 110]/255;
        case 'joint'
            model_colors(modelnum,:) = [38, 34, 98]/255;
        case 'markers'
            model_colors(modelnum,:) = [193, 25, 47]/255;
        end
    end

