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
    
    % handVel = muscVel*J'
    J_pm = muscVel_pm\handVel_pm;
    J_pm = J_pm';
    J_dl = muscVel_dl\handVel_dl;
    J_dl = J_dl';
    J_total = [muscVel_pm;muscVel_dl]\[handVel_pm;handVel_dl];
    J_total = J_total';
    
    % rearrange to only get x and y components
    J_pm = J_pm([3;1],:);
    J_dl = J_dl([3;1],:);
    J_total = J_total([3;1],:);
    
    % get change in Jacobian
    absJ_pm = sqrt(sum((J_pm).^2,1));
    angJ_pm = atan2d(J_pm(2,:),J_pm(1,:));
    absJ_dl = sqrt(sum((J_dl).^2,1));
    angJ_dl = atan2d(J_dl(2,:),J_dl(1,:));
    delJ = sum(sum((J_dl-J_pm).^2,1))
    shiftJ = mean(abs(angJ_dl-angJ_pm));

    % Get weighted change in Jacobian (weighted by eigenvalues from PCA)
    % This should also be the same as empirically calculating the difference
    %       in hand velocity between the two Jacobians given the same muscle
    %       velocities (assuming the muscle velocities are in keeping with
    %       PCA calculated above)
    % This is akin to a mean squared error between Jacobians
    muscVel_pca = get_vars([td_pm td_dl],{'opensim_muscVel_pca',1:5});
    N = length(muscVel_pca);
    delJ_weighted = (N-1)/N * sum((J_dl-J_pm).^2,1) * pca_info_vel.eigen(1:5)
    
    % Calculate empical difference in hand velocity between Jacobians
    handVel_pm_hat = muscVel_pca*J_pm';
    handVel_dl_hat = muscVel_pca*J_dl';
    del_handVel = mean(sum((handVel_dl_hat-handVel_pm_hat).^2,2),1)

    % Calculate empical difference in hand velocity between using unified Jacobian and actual hand velocity
    handVel_hat = muscVel_pca*J_total';
    handVel = [handVel_pm;handVel_dl];
    handVel = handVel(:,[3 1]);
    zeroidx = find(sum(handVel.^2,2)==0); % indices of dropped points most likely
    handVel_hat(zeroidx,:) = [];
    handVel(zeroidx,:) = [];
    del_handVel_frac = mean(sum((handVel_hat-handVel).^2,2)./sum(handVel.^2,2),1)

    % tuning analysis
    % for any neural tuning to muscle velocity (q_dot) that looks like 
    %   f = q_dot * \alpha + \epsilon_q (1)
    %   where f is firing rate and \alpha is muscle tuning parameters
    % there will be a corresponding endpoint tuning
    %   f = v_h * \beta + \epsilon_h (2)
    %   where f is firing rate, v_h is hand velocity, and \beta is endpoint tuning paramters
    % There exists a one-to-one transformation between q_dot and v_h, given a configuration q
    % (call this transformation v_h = J(q,q_dot), where J is linear in q_dot given q)
    % Assuming small excursions in q around a mean, J(q,q_dot) can be reduced to a linear transform
    % J, i.e. the average Jacobian matrix.
    %   v_h = q_dot * J' (3)
    % Thus, eq. 2 can be rewritten as
    %   f = q_dot * J' * \beta + \epsilon_h (4)
    % Combining (1) and (4)...
    %   q_dot * J' * \beta - \epsilon_h = q_dot * \alpha + \epsilon_q (5)
    %   q_dot * J' * \beta + \epsilon_c = q_dot * \alpha (6)
    % In a least squares sense...
    %   \alpha = J' * \beta
    % Are there \alpha's that would lead to identical or very similar \beta's?
    % Solve for \beta...
    %   \beta = ( (J * J')^-1 * J ) * \alpha
    % Taking the difference between two separate workspaces with possibly two different Jacobians...
    %   \Delta \beta = ( (J_dl * J_dl')^-1 * J_dl - (J_pm * J_pm')^-1 * J_pm ) * \alpha
    % Therefore, to find muscle tunings that cause little to no change in endpoint tuning, find the nullspace of
    %   (J_dl * J_dl')^-1 * J_dl - (J_pm * J_pm')^-1 * J_pm
    % Making \alpha p-dimensional, this nullspace is at least (p-2)-dimensional
    alpha_null = null( pinv(J_dl')-pinv(J_pm') )
    alpha_potent = orth( pinv(J_dl')'-pinv(J_pm')' )
