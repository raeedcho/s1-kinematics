%%
file_path = 'C:\Users\Raeed\NeuralData\Han_13B1\experiment_20151211_twoRW\RAW\';
file_prefix_all = 'Han_20141211_RW';

% Run processSpikesForSorting for the first time to combine spike data from
% all files with a name starting with file_prefix.
mergingStatus = processSpikesForSorting(file_path,file_prefix_all);

%% Check that the spike data has been successfully merged
while strcmp(mergingStatus,'merged spikes')
    % Now go to OfflineSorter and sort your spikes!
    disp(['Sort ''' file_prefix_all '-spikes.nev'' in OfflineSorter and save sorted file as '''...
        file_prefix_all '-spikes-s.nev'' then press any key to continue.'])
    pause
    % Run processSpiesForSorting again to separate sorted spikes into their
    % original files.
    mergingStatus = processSpikesForSorting(file_path,file_prefix_all);
    if strcmp(mergingStatus,'processed')
        % If everything went well, create bdfs for your files (you might
        % want to split them up by task.)
        bdf_DL = get_nev_mat_data([file_path file_prefix_all '_DL'],6);
        bdf_PM = get_nev_mat_data([file_path file_prefix_all '_PM'],6);
    end
end