function [alpha_potent,alpha_null] = get_alpha_potent(trial_data)
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
    PCAparams_vel = struct('signals',{{'opensim',find(contains(td_pm(1).opensim_names,'_muscVel'))}}, 'do_plot',false);
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
    
    % handVel = muscVel*J'
    J_pm = muscVel_pm\handVel_pm;
    J_pm = J_pm';
    J_dl = muscVel_dl\handVel_dl;
    J_dl = J_dl';
    
    % rearrange to only get x and y components
    J_pm = J_pm([3;1],:);
    J_dl = J_dl([3;1],:);

    % tuning analysis
    % for any neural tuning to muscle velocity (q_dot) that looks like 
    %   f = q_dot * \alpha + \epsilon_q (1)
    %   where f is firing rate and \alpha is muscle tuning parameters
    % there will be a corresponding endpoint tuning
    %   f = x_dot * \beta + \epsilon_h (2)
    %   where f is firing rate, x_dot is hand velocity, and \beta is endpoint tuning paramters
    % There exists a one-to-one transformation between q_dot and x_dot, given a configuration q
    % (call this transformation x_dot = J(q,q_dot), where J is linear in q_dot given q)
    % Assuming small excursions in q around a mean, J(q,q_dot) can be reduced to a linear transform
    % J, i.e. the average Jacobian matrix.
    %   x_dot = q_dot * J' (3)
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
    alpha_null = null( pinv(J_dl')-pinv(J_pm') );
    alpha_potent = orth( pinv(J_dl')'-pinv(J_pm')' );
