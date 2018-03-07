%% Estimate Jacobian for each epoch
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
    
    % get velocity PCA
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams_vel = struct('signals',{{'opensim',find(contains(td_pm(1).opensim_names,'_muscVel'))}}, 'do_plot',true);
    [~,pca_info_vel] = getPCA(cat(1,td_pm,td_dl),PCAparams_vel);
    td_pm = getPCA(td_pm,pca_info_vel);
    td_dl = getPCA(td_dl,pca_info_vel);
    % temporary hack to allow us to save into something useful
    for i=1:length(td_pm)
        td_pm(i).opensim_muscVel_pca = td_pm(i).opensim_pca;
        td_dl(i).opensim_muscVel_pca = td_dl(i).opensim_pca;
    end
    % get rid of superfluous PCA
    td_pm = rmfield(td_pm,'opensim_pca');
    td_dl = rmfield(td_dl,'opensim_pca');

    % get model of muscles to hand (fit Jacobian)
    % muscVel_pm = get_vars(td_pm,{'opensim',find(contains(td(1).opensim_names,'_muscVel'))});
    % muscVel_dl = get_vars(td_dl,{'opensim',find(contains(td(1).opensim_names,'_muscVel'))});
    muscVel_pm = get_vars(td_pm,{'opensim_muscVel_pca',1:5});
    muscVel_dl = get_vars(td_dl,{'opensim_muscVel_pca',1:5});
    handVel_pm = get_vars(td_pm,{'opensim',find(contains(td(1).opensim_names,'_handVel'))});
    handVel_dl = get_vars(td_dl,{'opensim',find(contains(td(1).opensim_names,'_handVel'))});
    handleVel_pm = get_vars(td_pm,{'vel',1:2})/100; % handleVel is in cm/s in TD
    handleVel_dl = get_vars(td_dl,{'vel',1:2})/100;
    
    % muscVel = handVel*J
    J_pm = handVel_pm\muscVel_pm;
    J_dl = handVel_dl\muscVel_dl;
    
    % rearrange to only get x and y components
    J_pm = J_pm([3;1],:);
    J_dl = J_dl([3;1],:);
    
    % get change in Jacobian
    absJ_pm = sqrt(sum((J_pm).^2,1));
    angJ_pm = atan2d(J_pm(2,:),J_pm(1,:));
    absJ_dl = sqrt(sum((J_dl).^2,1));
    angJ_dl = atan2d(J_dl(2,:),J_dl(1,:));
    delJ = sum(sqrt(sum((J_dl-J_pm).^2,1)))
    shiftJ = mean(abs(angJ_dl-angJ_pm))
