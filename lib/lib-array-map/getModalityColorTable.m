function modality_color_table = getModalityColorTable()
% Function to get a table of all colors corresponding to a sensory modality

% modality colors from color brewer (diverging)
modality_colors = [...
    252,141,89;...
    153,213,148]/255;
modalities = {'cutaneous';'deep'};
modality_color_table = table(modalities,modality_colors,'VariableNames',{'modality','modality_color'});
