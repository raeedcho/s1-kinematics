%% Script to put together preliminary data for Kramer PD two workspace

ul = unit_list(bdf1,1);

%% loop over units
b_center = zeros(6,length(ul));
b_side = b_center;
PD_angle_diff = zeros(length(ul),1);

%%
tic;
for i = 1:length(ul)
    et = toc;
    fprintf(1, 'ET: %f (%d of %d)\n', et, i, length(ul));
    
    %% Set up model inputs
    ts = 50;

    vt_1 = bdf1.vel(:,1);
    vt_2 = bdf2.vel(:,1);
    vt_3 = bdf3.vel(:,1);
    
    t_1 = vt_1(1):ts/1000:vt_1(end);
    t_2 = vt_2(1):ts/1000:vt_2(end);
    t_3 = vt_3(1):ts/1000:vt_3(end);
    
    spike_times_1 = get_unit(bdf1,ul(i,1),ul(i,2));
    spike_times_2 = get_unit(bdf2,ul(i,1),ul(i,2));
    spike_times_3 = get_unit(bdf3,ul(i,1),ul(i,2));
    
    spike_times_1 = spike_times_1(spike_times_1>t_1(1) & spike_times_1<t_1(end));
    spike_times_2 = spike_times_2(spike_times_2>t_2(1) & spike_times_2<t_2(end));
    spike_times_3 = spike_times_3(spike_times_3>t_3(1) & spike_times_3<t_3(end));
    
    s_1 = train2bins(spike_times_1, t_1);
    s_2 = train2bins(spike_times_2, t_2);
    s_3 = train2bins(spike_times_3, t_3);
    
    % bdf1 and bdf2 are both center workspace
    s_center = [s_1 s_2];
    s_side = s_3;

%     % set default value of num_samp
%     if(isNaN(num_samp))
%         num_samp = length(t);
%     end
    
    glmx1 = interp1(bdf1.pos(:,1), bdf1.pos(:,2:3), t_1);
    glmv1 = interp1(bdf1.vel(:,1), bdf1.vel(:,2:3), t_1);
    glmx2 = interp1(bdf2.pos(:,1), bdf2.pos(:,2:3), t_2);
    glmv2 = interp1(bdf2.vel(:,1), bdf2.vel(:,2:3), t_2);
    glmx3 = interp1(bdf3.pos(:,1), bdf3.pos(:,2:3), t_3);
    glmv3 = interp1(bdf3.vel(:,1), bdf3.vel(:,2:3), t_3);

    % bdf1 and bdf2 are both center workspace
    glm_input_center = [glmx1 glmv1 sqrt(glmv1(:,1).^2 + glmv1(:,2).^2);...
                        glmx2 glmv2 sqrt(glmv2(:,1).^2 + glmv2(:,2).^2)];
    glm_input_side = [glmx3 glmv3 sqrt(glmv3(:,1).^2 + glmv3(:,2).^2)];
    
%     %% Bootstrap
%     % bootstrap
%     b_mat = zeros(size(glm_input,2)+1,1,reps);
%     for bootCt=1:reps
%         % grab test set indices
%         idx = uint32(1+(length(glmx)-1)*rand(num_samp,1));
% 
%         b = glmfit(glm_input(idx,:),s(idx),'poisson');
%         moddepth_boot(i,bootCt) = norm([b(4) b(5)]);
%         pds_boot(i,bootCt) = atan2(b(5),b(4));
%         
%         b_mat(:,bootCt) = b;
%     end
%     avg_b = mean(b_mat,2);

    % find indices to keep
    center_idx = [2250:6200 6600:7900 8000:14720];
    side_idx = 1:1110;
    
    glm_input_center = glm_input_center(center_idx,:);
    glm_input_side = glm_input_side(side_idx,:);
    s_center = s_center(center_idx);
    s_side = s_side(side_idx);
    
%     center_idx = [225:620 660:790 800:1472];
%     side_idx = 1:111;

    %% GLM fit both workspaces
    for repct = 1:1000
        center_idx_rand = uint32(1+(length(center_idx)-1)*rand(length(center_idx),1));
        side_idx_rand = uint32(1+(length(side_idx)-1)*rand(length(side_idx),1));
        b_bootstrap_center(:,repct) = glmfit(glm_input_center(center_idx_rand,:),s_center(center_idx_rand),'poisson');
        b_bootstrap_side(:,repct) = glmfit(glm_input_side(side_idx_rand,:),s_side(side_idx_rand),'poisson');
        PD_boot_center(repct) = atan2d(b_bootstrap_center(5,repct),b_bootstrap_center(4,repct));
        PD_boot_side(repct) = atan2d(b_bootstrap_side(5,repct),b_bootstrap_side(4,repct));
        
        % Do full regression
        full_glm_input = [glm_input_center zeros(length(glm_input_center),size(glm_input_side,2)+1); glm_input_side ones(length(glm_input_side),1) glm_input_side];
        s_full = [s_center s_side]';
        full_idx_rand = uint32(1+(length(full_glm_input)-1)*rand(length(full_glm_input),1));
        b_bootstrap_full = GeneralizedLinearModel.fit(full_glm_input(full_idx_rand,:),s_full(full_idx_rand),'Distribution','poisson');
        b_full_est(:,repct) = b_bootstrap_full.Coefficients.Estimate;
    end
    
    b_full_mean(:,i) = mean(b_full_est,2);
    b_full_sort = sort(b_full_est,2,'ascend');
    b_full_min(:,i) = b_full_sort(:,10);
    b_full_max(:,i) = b_full_sort(:,990);
    
    b_center(:,i) = mean(b_bootstrap_center,2);
    b_side(:,i) = mean(b_bootstrap_side,2);
    
    b_center_deet{i} = GeneralizedLinearModel.fit(glm_input_center,s_center,'Distribution','poisson');
    b_side_deet{i} = GeneralizedLinearModel.fit(glm_input_side,s_side,'Distribution', 'poisson');
    
    
    % calculate PD change
    PD_angle_diff(i) = acosd((b_center(2:6,i)'*b_side(2:6,i))/(norm(b_center(2:6,i))*norm(b_side(2:6,i))));
%     PD_angle_diff_simple(i) = acosd((b_center(4:5,i)'*b_side(4:5,i))/(norm(b_center(4:5,i))*norm(b_side(4:5,i))));
    PD_simple_center(i) = atan2d(b_center(5,i),b_center(4,i));
    PD_simple_side(i) = atan2d(b_side(5,i),b_side(4,i));
    PD_angle_diff_simple(i) = PD_simple_side(i)-PD_simple_center(i);
    
    % find CI
    PD_sort_center = sort(mod(PD_boot_center-PD_simple_center(i)+180,360)-180);
    PD_sort_side = sort(mod(PD_boot_side-PD_simple_side(i)+180,360)-180);
    
    CI_center(1,i) = PD_sort_center(25)+PD_simple_center(i);
    CI_center(2,i) = PD_sort_center(975)+PD_simple_center(i);
    CI_side(1,i) = PD_sort_side(25)+PD_simple_side(i);
    CI_side(2,i) = PD_sort_side(97)+PD_simple_side(i);
    
    CI_range_center(i) = CI_center(2,i)-CI_center(1,i);
    CI_range_side(i) = CI_side(2,i)-CI_side(1,i);
    moddepth_center(i) = norm(b_center(4:5,i));
    moddepth_side(i) = norm(b_side(4:5,i));
end

%%

sig_change = ((b_full_min>0 & b_full_mean>0) | (b_full_max<0 & b_full_mean<0));
good_units = CI_range_center<=45 & CI_range_side<=45;
sig_change_good = sig_change(:,good_units);
sig_change_good = sig_change(7:12,good_units);
sig_collapse = sum(sig_change);
sig_collapse = sum(sig_change_good);
sig_bad = sig_change(:,~good_units);
sig_count = sum(sig_change_good,2);

PD_angle_diff_simple = mod(PD_angle_diff_simple+180,360)-180;

% figure;hist(PD_angle_diff,20);
figure; hist(PD_angle_diff_simple(good_units),15);
title 'Histogram of PD changes across workspaces'
xlabel '\Delta PD'
ylabel 'Number of neurons'

figure;hist(PD_simple_center(good_units))
title 'Histogram of PD in blue workspace'
xlabel 'PD'
ylabel 'Number of neurons'
figure;hist(PD_simple_side(good_units))
title 'Histogram of PD in red workspace'
xlabel 'PD'
ylabel 'Number of neurons'


%%
figure; 
h1 = plot(glm_input_center(1:100,1),glm_input_center(1:100,2),'-b',glm_input_side(1:100,1),glm_input_side(1:100,2),'-r');
set(h1,'LineWidth',3)

%%
% figures of PDs (37)
figure;
h4 = polar(0,0.05,'.');
hold on
h2 = polar([PD_simple_center(41) PD_simple_center(41)],[0 moddepth_center(41)],'-b');
h3 = polar(PD_simple_center(41),moddepth_center(41),'ob');
set(h4,'MarkerSize',1)
set(h2,'LineWidth',3)
set(h3,'LineWidth',3,'MarkerSize',10)
h2 = polar([PD_simple_side(41) PD_simple_side(41)],[0 moddepth_side(41)],'-r');
h3 = polar(PD_simple_side(41),moddepth_side(41),'or');
set(h4,'MarkerSize',1)
set(h2,'LineWidth',3)
set(h3,'LineWidth',3,'MarkerSize',10)

