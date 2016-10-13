%% Script to put together preliminary data for Kramer PD two workspace

ul = unit_list(bdf1,1);

%% loop over units
tic;
for i = 1:length(ul)
    et = toc;
    fprintf(1, 'ET: %f (%d of %d)\n', et, i, length(ul));
    
    %% Set up model inputs
    ts = 50;

    vt_1 = bdf1.vel(:,1);
    vt_2 = bdf2.vel(:,1);
    vt_3 = bdf3.vel(:,1);
    vt_4 = bdf4.vel(:,1);
    
    t_1 = vt_1(1):ts/1000:vt_1(end);
    t_2 = vt_2(1):ts/1000:vt_2(end);
    t_3 = vt_3(1):ts/1000:vt_3(end);
    t_4 = vt_4(1):ts/1000:vt_4(end);
    
    spike_times_1 = get_unit(bdf1,ul(i,1),ul(i,2));
    spike_times_2 = get_unit(bdf2,ul(i,1),ul(i,2));
    spike_times_3 = get_unit(bdf3,ul(i,1),ul(i,2));
    spike_times_4 = get_unit(bdf4,ul(i,1),ul(i,2));
    
    spike_times_1 = spike_times_1(spike_times_1>t_1(1) & spike_times_1<t_1(end));
    spike_times_2 = spike_times_2(spike_times_2>t_2(1) & spike_times_2<t_2(end));
    spike_times_3 = spike_times_3(spike_times_3>t_3(1) & spike_times_3<t_3(end));
    spike_times_4 = spike_times_4(spike_times_4>t_4(1) & spike_times_4<t_4(end));
    
    s_1 = train2bins(spike_times_1, t_1);
    s_2 = train2bins(spike_times_2, t_2);
    s_3 = train2bins(spike_times_3, t_3);
    s_4 = train2bins(spike_times_4, t_4);
    
    % bdf1 and bdf2 are both center workspace
    s = [s_1 s_2 s_3 s_4];

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
    glmx4 = interp1(bdf4.pos(:,1), bdf4.pos(:,2:3), t_4);
    glmv4 = interp1(bdf4.vel(:,1), bdf4.vel(:,2:3), t_4);

    % bdf1 and bdf2 are both center workspace
    glm_input = [glmx1 glmv1 sqrt(glmv1(:,1).^2 + glmv1(:,2).^2);...
                 glmx2 glmv2 sqrt(glmv2(:,1).^2 + glmv2(:,2).^2);...
                 glmx3 glmv3 sqrt(glmv3(:,1).^2 + glmv3(:,2).^2);...
                 glmx4 glmv4 sqrt(glmv4(:,1).^2 + glmv4(:,2).^2)];

    % find indices to keep
    keep_idx = [1:14150 14250:19100 19400:22700 23000:28650 28950:30500 30700:32600 330500:35700 36200:37100 37600:38200 38500:40900 41000:41800 42100:42900 43100:43900 44100:47400 47600:48100 48400:49100 49400:52700 52900:61800 62200:64600 64750:66200 66400:66650 66850:68950 69350:71950 72100:74200 75100:75400 75900:77300 77700:78750 79700:80200 80500:81850 82400:82900 83200:83900 84900:85900 87100:88450 89000:91900 93300:94600 94900:97500 97750:102400 103600:105300 105700:108700 110200:113500 114200:121700];
    
    glm_input = glm_input(keep_idx,:);
    s = s(keep_idx);
    
%     center_idx = [225:620 660:790 800:1472];
%     side_idx = 1:111;

    %% GLM fit both workspaces
    for repct = 1:1000
        idx_rand = uint32(1+(length(keep_idx)-1)*rand(length(keep_idx),1));
        b_bootstrap(:,repct) = glmfit(glm_input(idx_rand,:),s(idx_rand),'poisson');
        PD_boot(repct) = atan2d(b_bootstrap(5,repct),b_bootstrap(4,repct));
    end
    
    b_est(:,i) = mean(b_bootstrap,2);
    
    % find CI
    PD_simple(i) = atan2d(b_est(5,i),b_est(4,i));
    PD_sort_boot = sort(mod(PD_boot-PD_simple(i)+180,360)-180);
    
    CI_boot(1,i) = PD_sort_boot(25)+PD_simple(i);
    CI_boot(2,i) = PD_sort_boot(975)+PD_simple(i);
    
    CI_range(i) = CI_boot(2,i)-CI_boot(1,i);
    moddepth(i) = norm(b_est(4:5,i));
end
