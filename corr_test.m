%% Preprocess trial_data
    % prep trial data by getting only rewards and trimming to only movements
    [~,td] = getTDidx(trial_data,'result','R');
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % bin data at 50ms
    td = binTD(td,5);
    
    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);

    % Bootstrap to find pairwise correlations
    num_boots = 1000;
    trial_idx = randi(minsize,minsize,num_boots);
    rho_pm = zeros(size(td(1).S1_spikes,2),size(td(1).S1_spikes,2),num_boots);
    rho_dl = zeros(size(td(1).S1_spikes,2),size(td(1).S1_spikes,2),num_boots);
    tic
    for bootctr = 1:num_boots
        [temp_rho_pm,sort_pm] = pairwiseCorr(td_pm(trial_idx(:,bootctr)),struct('signals',{{'S1_spikes'}},'cluster_order',true));
        [temp_rho_dl,sort_dl] = pairwiseCorr(td_dl(trial_idx(:,bootctr)),struct('signals',{{'S1_spikes'}},'cluster_order',true));
        rho_pm(:,:,bootctr) = temp_rho_pm;
        rho_dl(:,:,bootctr) = temp_rho_dl;
        disp(['Bootstrap sample ' num2str(bootctr) ', done at ' num2str(toc) 's'])
    end
    rho_diff = rho_dl-rho_pm;
    
    mean_rho_pm = mean(rho_pm,3);
    mean_rho_dl = mean(rho_dl,3);
    mean_rho_diff = mean(rho_diff,3);
    
    CI_rho_pm = prctile(rho_pm,[2.5 97.5],3);
    CI_rho_dl = prctile(rho_dl,[2.5 97.5],3);
    CI_rho_diff = prctile(rho_diff,[2.5 97.5],3);
    
    rho_diff_zeros = CI_rho_diff(:,:,1)<0 & CI_rho_diff(:,:,2)>0;
    rho_diff_contrast = mean_rho_diff;
    rho_diff_contrast(rho_diff_zeros) = 0;
    
    clim = [min(min([mean_rho_dl mean_rho_pm])) max(max([mean_rho_dl mean_rho_pm]))];

    % display
    figure
    subplot(1,2,1)
    imagesc(mean_rho_pm(sort_pm,sort_pm),clim)
    subplot(1,2,2)
    imagesc(mean_rho_dl(sort_pm,sort_pm),clim)
    figure
    subplot(1,2,1)
    imagesc(mean_rho_diff(sort_pm,sort_pm),clim)
    subplot(1,2,2)
    imagesc(rho_diff_contrast(sort_pm,sort_pm),clim)
    
    % display
    figure
    subplot(1,2,1)
    mesh(mean_rho_pm(sort_pm,sort_pm))
    set(gca,'zlim',clim)
    subplot(1,2,2)
    mesh(mean_rho_dl(sort_pm,sort_pm))
    set(gca,'zlim',clim)
    figure
    mesh(mean_rho_diff(sort_pm,sort_pm))
    set(gca,'zlim',clim)
    
    % clean up
    clearvars td td_pm td_dl minsize *_idx
