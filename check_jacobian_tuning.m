%% Get alphas from glm_info
    % Quick analysis
    results_quick = analyzeTRT(trial_data,struct('num_boots',1));

    % get alpha
    alpha = results_quick.glm_info{3}.b(7:11,:);
    alpha_mag = sqrt(sum(alpha.^2,1));

    % get alpha potent space
    alpha_potent = get_alpha_potent(trial_data);

    % Get projection
    alpha_proj = alpha_potent'*alpha;
    alpha_proj_mag = sqrt(sum(alpha_proj.^2,1));
    potent_frac = alpha_proj_mag./alpha_mag;

    % histogram
    figure
    hist(potent_frac)
    set(gca,'box','off','tickdir','out','xlim',[0,1])
