%% Just plot modulation for one channel
% set channel
chan = 55;

% get unit list
ul = unit_list(bdf1);

% get only units of interest
ul = ul(ul(:,1)==chan,:);

%% Find instantaneous firing rates of units
for i = 1:length(ul)
    ts = 50;

    vt = bdf1.vel(:,1);
    t = vt(1):ts/1000:vt(end);
    spike_times = get_unit(bdf1,ul(i,1),ul(i,2));
    spike_times = spike_times(spike_times>t(1) & spike_times<t(end));
    s = train2bins(spike_times, t);
    
    b = glm_kin(bdf1,ul(i,1),ul(i,2),0,'posvel');
    PD_vec = [b(4);b(5)];
    unit_PD_vec = PD_vec/norm(PD_vec);
    
    % convolve with gaussian
    gauss_filt = gausswin(6);
    gauss_filt = gauss_filt/sum(gauss_filt);
    f_est(:,i) = transpose(conv(s,gauss_filt,'same'));
    speed_temp = bdf1.vel(:,2:3)*unit_PD_vec;
    speed(:,i) = interp1(vt,speed_temp,t);
end

%% plot
subplot(2,1,1)
plot(t(1:600),20*f_est(1:600,:))
title 'Estimated Firing Rate'
ylabel 'Estimated firing rate (Hz)'
subplot(2,1,2)
plot(t(1:600),speed(1:600,:))
title 'Velocity along preferred direction'
legend('Channel 55, Unit 1','Channel 55, Unit 2', 'Channel 55, Unit 3')
xlabel 'Time (s)'
ylabel 'Velocity (cm/s)'