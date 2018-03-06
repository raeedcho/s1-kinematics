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

    % Get pairwise correlations
    % [rho_pm,sort_pm] = pairwiseCorr(td_pm,struct('signals',{{'S1_spikes'}},'cluster_order',true));
    [rho_pm,sort_pm] = pairwiseCorr(td_pm,struct('signals',{{'S1_spikes'}}));
    [rho_dl,sort_dl] = pairwiseCorr(td_pm,struct('signals',{{'S1_spikes'}}));
    rho_dl  = rho_dl(sort_pm,sort_pm);
    clim = [min(min([rho_dl rho_pm])) max(max([rho_dl rho_pm]))];

    % display
    figure
    subplot(1,2,1)
    imagesc(rho_pm,clim)
    subplot(1,2,2)
    imagesc(rho_dl,clim)
    
    % clean up
    clearvars td td_pm td_dl minsize *_idx
