%% Prep TD
    % prep trial data by getting only rewards and trimming to only movements
    % first process marker data
    td = trial_data;
    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signal','markers','alias','marker_vel'));
    [~,td] = getTDidx(td,'result','R');
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % for bumps
    % [~,td] = getTDidx(trial_data,'result','R');
    % td = td(~isnan(cat(1,td.idx_bumpTime)));
    % td = trimTD(td,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % bin data at 50ms
    td = binTD(td,5);
    % add in spherical coordinates
    td = addSphereHand2TD(td);
    % add in cylindrical coordinates
    td = addCylHand2TD(td);
    % add firing rates rather than spike counts
    td = addFiringRates(td,struct('array','S1'));

    %% Do PCA on muscle space
    % do PCA on muscles, training on only the training set
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len'))}}, 'do_plot',true);
    [td,~] = getPCA(td,PCAparams);
    % temporary hack to allow us to do PCA on velocity too
    for i=1:length(td)
        td(i).opensim_len_pca = td(i).opensim_pca;
    end
    % get rid of superfluous PCA
    td = rmfield(td,'opensim_pca');
    % get velocity PCA
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel'))}}, 'do_plot',true);
    [td,~] = getPCA(td,PCAparams_vel);
    % temporary hack to allow us to save into something useful
    for i=1:length(td)
        td(i).opensim_muscVel_pca = td(i).opensim_pca;
    end
    % get rid of superfluous PCA
    td = rmfield(td,'opensim_pca');

    % Get PCA for neural space
    % PCAparams = struct('signals',{{'S1_spikes'}}, 'do_plot',true,'pca_recenter_for_proj',true,'sqrt_transform',true);
    % [td,~] = getPCA(td,PCAparams);

    % Get PCA for marker space
    td = getPCA(td,struct('signals','markers'));
    td = getPCA(td,struct('signals','marker_vel'));

    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);

    % recombine for later...
    td = [td_pm td_dl];

%% Create model from marker pca to muscle pca
    marker_frame = get_vars(td,{'markers_pca',1:5;'marker_vel_pca',1:5});
    muscle_frame = get_vars(td,{'opensim_len_pca',1:5;'opensim_muscVel_pca',1:5});

    model = marker_frame\muscle_frame;
    muscle_pred = marker_frame*model;
    muscle_resid = muscle_frame-muscle_pred;

%% test model linearity by plotting predictions and residuals
    for i = 1:size(muscle_frame,2)
        figure
        subplot(1,2,1)
        scatter(muscle_frame(:,i),muscle_pred(:,i),[],'k','filled')
        hold on
        plot([min(muscle_frame(:,i)) max(muscle_frame(:,i))],[min(muscle_frame(:,i)) max(muscle_frame(:,i))],'--k','linewidth',2)
        axis square
        set(gca,'box','off','tickdir','out')
        xlabel(sprintf('Muscle frame PC %d',i))
        ylabel(sprintf('Predicted muscle frame PC %d',i))
        title(sprintf('Prediction plot for PC %d',i))

        subplot(1,2,2)
        scatter(muscle_frame(:,i),muscle_resid(:,i),[],'k','filled')
        hold on
        plot([min(muscle_frame(:,i)) max(muscle_frame(:,i))],[0 0],'--k','linewidth',2)
        axis square
        set(gca,'box','off','tickdir','out')
        xlabel(sprintf('Muscle frame PC %d',i))
        ylabel(sprintf('Predicted muscle frame PC %d',i))
        title(sprintf('Residual plot for PC %d',i))
    end

    figure
    scatter3(muscle_frame(:,1),muscle_frame(:,2),muscle_frame(:,3),[],'k','filled')
    axis equal
    figure
    scatter3(muscle_pred(:,1),muscle_pred(:,2),muscle_pred(:,3),[],'k','filled')
    axis equal

    err = sqrt(sum(muscle_resid(:,1:2).^2,2));
    figure
    scatter3(muscle_frame(:,1),muscle_frame(:,2),err,[],'k','filled')
    axis equal

%% neural network?
    net = feedforwardnet([10 10]);
    net = train(net,marker_frame',muscle_frame');

    net_pred = net(marker_frame');
    net_pred = net_pred';
    net_resid = muscle_frame - net_pred;

    for i = 1:size(muscle_frame,2)
        figure
        subplot(1,2,1)
        scatter(muscle_frame(:,i),net_pred(:,i),[],'k','filled')
        hold on
        plot([min(muscle_frame(:,i)) max(muscle_frame(:,i))],[min(muscle_frame(:,i)) max(muscle_frame(:,i))],'--k','linewidth',2)
        axis square
        set(gca,'box','off','tickdir','out')
        xlabel(sprintf('Muscle frame PC %d',i))
        ylabel(sprintf('Predicted muscle frame PC %d',i))
        title(sprintf('Prediction plot for PC %d',i))

        subplot(1,2,2)
        scatter(muscle_frame(:,i),net_resid(:,i),[],'k','filled')
        hold on
        plot([min(muscle_frame(:,i)) max(muscle_frame(:,i))],[0 0],'--k','linewidth',2)
        axis square
        set(gca,'box','off','tickdir','out')
        xlabel(sprintf('Muscle frame PC %d',i))
        ylabel(sprintf('Predicted muscle frame PC %d',i))
        title(sprintf('Residual plot for PC %d',i))
    end
