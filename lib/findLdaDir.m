function true_lda_coeff = findLdaDir(lda_table)
    % This function finds the direction of the LDA axis in the full neural space

    % set up key columns in final table
    keycols = strcmpi(lda_table.Properties.VariableDescriptions,'meta');
    key_table = lda_table(:,keycols);

    % for each entry, multiply lda vectory by pca weights to get lda vector in neural space
    coeff_mat = zeros(height(key_table),length(lda_table(1,:).S1_FR_pca_coeff{1}));
    for entrynum = 1:height(key_table)
        entry_coeff = lda_table(entrynum,:).S1_FR_pca_coeff{1}*lda_table(entrynum,:).S1_FR_lda_coeff(2:end)';
        coeff_mat(entrynum,:) = entry_coeff;
    end

    % package into NeuronTable
    true_lda_coeff = horzcat(key_table,table(coeff_mat,'VariableNames',{'S1_FR_lda_full_coeff'}));
end
