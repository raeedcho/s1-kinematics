function [figure_handles, output_data]=rotate_tuning_curves(folder,options)
% Find scaling and rotation between two tuning curves in PM and DL
% workspaces

%% Get data to fit
% % prep bdf
% options.bdf_PM.meta.task = 'RW';
% options.bdf_DL.meta.task = 'RW';
% 
% %add firing rate to the units fields of the bdf
% opts.binsize=0.05;
% opts.offset=-.015;
% opts.do_trial_table=1;
% opts.do_firing_rate=1;
% options.bdf_PM=postprocess_bdf(options.bdf_PM,opts);
% options.bdf_DL=postprocess_bdf(options.bdf_DL,opts);
% 
% optionstruct.compute_pos_pds=0;
% optionstruct.compute_vel_pds=1;
% optionstruct.compute_acc_pds=0;
% optionstruct.compute_force_pds=0;
% optionstruct.compute_dfdt_pds=0;
% optionstruct.compute_dfdtdt_pds=0;
% if(isfield(options,'which_units'))
%     which_units = options.which_units;
% else
%     for i=1:length(options.bdf_PM.units)
%         temp(i)=options.bdf_PM.units(i).id(2)~=0 && options.bdf_PM.units(i).id(2)~=255;
%     end
%     ulist=1:length(options.bdf_PM.units);
%     which_units=ulist(temp);
% end
% optionstruct.data_offset=-.015;%negative shift shifts the kinetic data later to match neural data caused at the latency specified by the offset
% 
% optionstruct_DL.compute_pos_pds=0;
% optionstruct_DL.compute_vel_pds=1;
% optionstruct_DL.compute_acc_pds=0;
% optionstruct_DL.compute_force_pds=0;
% optionstruct_DL.compute_dfdt_pds=0;
% optionstruct_DL.compute_dfdtdt_pds=0;
% if(isfield(options,'which_units'))
%     which_units = options.which_units;
% else
%     for i=1:length(options.bdf_DL.units)
%         temp(i)=options.bdf_DL.units(i).id(2)~=0 && options.bdf_DL.units(i).id(2)~=255;
%     end
%     ulist=1:length(options.bdf_DL.units);
%     which_units_DL=ulist(temp);
% end
% optionstruct_DL.data_offset=-.015;%negative shift shifts the kinetic data later to match neural data caused at the latency specified by the offset


% set up PM tuning curve function
options_PM = options;
options_PM.plot_curves=0;
if isfield(options,'bdf_PM')
    options_PM.bdf = options.bdf_PM;
    options_PM = rmfield(options_PM,'bdf_PM');
end
    
% behaviors_PM = parse_for_tuning(options.bdf_PM,'continuous','opts',optionstruct,'units',which_units);
% options_PM.behaviors = behaviors_PM;

[~,tuning_PM] = get_tuning_curves(folder,options_PM);

% set up DL tuning curve function
options_DL = options;
options_DL.plot_curves=0;
if isfield(options,'bdf_DL')
    options_DL.bdf = options.bdf_DL;
    options_DL = rmfield(options_DL,'bdf_DL');
end
% behaviors_DL = parse_for_tuning(options.bdf_DL,'continuous','opts',optionstruct_DL,'units',which_units_DL);
% options_DL.behaviors = behaviors_DL;

[~,tuning_DL] = get_tuning_curves(folder,options_DL);

% get relevant data
binned_spd_PM = tuning_PM.binned_spd;
binned_spd_DL = tuning_DL.binned_spd;
FR_PM = tuning_PM.binned_FR;
FR_DL = tuning_DL.binned_FR;
angs_PM = repmat(tuning_PM.bins,1,size(FR_PM,2));
angs_DL = repmat(tuning_DL.bins,1,size(FR_DL,2));

%% Convert to complex polar representation
polar_PM_curve = FR_PM.*exp(1i*angs_PM);
polar_DL_curve = FR_DL.*exp(1i*angs_DL);

%% Fit data
for i = 1:size(FR_PM,2)
%     tbl = table(polar_PM_curve(:,i),polar_DL_curve(:,i),'VariableNames',{'PM_curve','DL_curve'});
%     lm = fitlm(tbl,'DL_curve ~ PM_curve - 1');
%     complex_scale_factor(1,i) = lm.Coefficients.Estimate;
%     
%     complex_scale_manual(1,i) = (polar_PM_curve(:,i)'*polar_PM_curve(:,i))\polar_PM_curve(:,i)'*(polar_DL_curve(:,i));

    % Fit with optimization
    real_imag_scale = fminsearch(@(x) (find_ms_curve_dist(polar_DL_curve(:,i),(x(1)+1i*x(2))*polar_PM_curve(:,i)))^2,rand(2,1));
    complex_scale_factor(1,i) = real_imag_scale(1)+1i*real_imag_scale(2);
    
    disp(['Done with ' num2str(i)])
end

scale_factor = abs(complex_scale_factor);
rot_factor = angle(complex_scale_factor);

%% Plot fits
figure_handles = [];
unit_ids = tuning_PM.unit_ids;
for i = 1:size(FR_PM,2)
    fig = figure('name',['channel_' num2str(unit_ids(i,1)) '_unit_' num2str(unit_ids(i,2)) '_tuning_plot']);
    figure_handles = [figure_handles fig];
    
    % polar tuning curve
    subplot(211)
    max_rad = max(abs([polar_PM_curve(:,i);polar_DL_curve(:,i)]));
    h=polar(0,max_rad);
    set(h,'color','w')
    hold on
    h=polar(angle(repmat(polar_PM_curve(:,i),2,1)),abs(repmat(polar_PM_curve(:,i),2,1)));
    set(h,'linewidth',2,'color',[0.6 0.5 0.7])
    hold on
    h=polar(angle(repmat(polar_DL_curve(:,i),2,1)),abs(repmat(polar_DL_curve(:,i),2,1)));
    set(h,'linewidth',2,'color',[1 0 0])
    h=polar(angle(repmat(polar_PM_curve(:,i)*complex_scale_factor(i),2,1)),abs(repmat(polar_PM_curve(:,i)*complex_scale_factor(i),2,1)));
    set(h,'linewidth',2,'color',[0 1 0])
    title 'Wrapped tuning curves'
    
    % flat tuning curve
    subplot(212)
    [rays_PM,mags_PM] = get_full_curve(polar_PM_curve(:,i));
    [rays_DL,mags_DL] = get_full_curve(polar_DL_curve(:,i));
    [rays_fit,mags_fit] = get_full_curve(polar_PM_curve(:,i)*complex_scale_factor(i));
    h=plot(180/pi*rays_PM,mags_PM);
    set(h,'linewidth',2,'color',[0.6 0.5 0.7])
    hold on
    h=plot(180/pi*rays_DL,mags_DL);
    set(h,'linewidth',2,'color',[1 0 0])
    h=plot(180/pi*rays_fit,mags_fit);
    set(h,'linewidth',2,'color',[0 1 0])
    
    set(gca,'xlim',[-180,180],'xtick',[-180 -90 0 90 180],'tickdir','out','box','off');
    xlabel 'Movement direction (deg)'
    ylabel 'Average spikes per 50 ms time bin'
    legend('PM curve','DL curve','Rotated/Scaled PM curve')
    legend('boxoff')
    title 'Unwrapped tuning curves'
    
%     saveas(fig,['channel_' num2str(unit_ids(i,1)) '_unit_' num2str(unit_ids(i,2)) '_tuning_plot' '.png'])
end

%% Plot binned speed
subplot(211)
max_rad = max(abs([binned_spd_PM;binned_spd_DL]));
h=polar(0,max_rad);
set(h,'color','w')
hold on
h=polar(repmat(tuning_PM.bins,2,1),repmat(binned_spd_PM,2,1));
set(h,'linewidth',2,'color',[0.6 0.5 0.7])
hold on
h=polar(repmat(tuning_DL.bins,2,1),repmat(binned_spd_DL,2,1));
set(h,'linewidth',2,'color',[1 0 0])
title 'Wrapped speed curves'

% flat tuning curve
subplot(212)
h=plot(tuning_PM.bins,binned_spd_PM);
set(h,'linewidth',2,'color',[0.6 0.5 0.7])
hold on
h=plot(tuning_DL.bins,binned_spd_DL);
set(h,'linewidth',2,'color',[1 0 0])

%% Plot summary
array_break_idx = find(unit_ids(:,1)>options.array_break,1,'first');

frac_moddepth_PM = tuning_PM.frac_moddepth;
frac_moddepth_DL = tuning_DL.frac_moddepth;
minfrac_moddepth = min([frac_moddepth_PM;frac_moddepth_DL]);

fig = figure('name','total_transform_plot');
figure_handles = [figure_handles fig];
max_rad = max(minfrac_moddepth);
h=polar(0,max_rad);
set(h,'color','w')
hold on

h=polar(repmat(rot_factor(1:array_break_idx-1),2,1),[zeros(size(minfrac_moddepth(1:array_break_idx-1)));minfrac_moddepth(1:array_break_idx-1)]);
set(h,'linewidth',2,'color',[0 0 1])
h=polar(repmat(rot_factor(array_break_idx:end),2,1),[zeros(size(minfrac_moddepth(array_break_idx:end)));minfrac_moddepth(array_break_idx:end)]);
set(h,'linewidth',2,'color',[0 1 0])
title({'Tuning curve rotations per area (Blue=3a, Green=2),', 'scaled by fractional modulation depth'})

saveas(fig,'total_transform_plot.png')

output_data.scale_factor = scale_factor;
output_data.rot_factor = rot_factor;

%% hypothesis testing
 % Check if change over each array is significantly different
% from zero
[hyp1,p1] = ttest(rot_factor(1:array_break_idx-1));
[hyp2,p2] = ttest(rot_factor(array_break_idx:end));
mean1 = mean(rot_factor(1:array_break_idx-1));
mean2 = mean(rot_factor(array_break_idx:end));
output_data.array_change_signif = [hyp1;hyp2];
output_data.array_change_pval = [p1;p2];
output_data.array_change_mean = [mean1;mean2];

% Test if changes over arrays are different from each other
[hyp,p] = ttest2(rot_factor(1:array_break_idx-1),rot_factor(array_break_idx:end),'Vartype','unequal');
output_data.between_array_signif = hyp;
output_data.between_array_pval = p;
output_data.between_array_mean = mean1-mean2;

output_data.num_units = [array_break_idx-1; length(rot_factor)-array_break_idx+1];

end

function [dist] = find_ms_curve_dist(curve1,curve2)
    %% add up distances along many rays
    % pad curves with wrap-around term for full curve
    curve1_wrap = [curve1(end);curve1];
    curve2_wrap = [curve2(end);curve2];
    
    % find full curves and remove duplicate point
    interpolater = linspace(1,length(curve1_wrap),length(curve1_wrap)*500+1)';
    interpolater = interpolater(2:end);
    curve1_full = interp1(curve1_wrap,interpolater);
    curve2_full = interp1(curve2_wrap,interpolater);
    
    % get sorted curves in angles and magnitudes
    curve1_sort = sortrows([angle(curve1_full) abs(curve1_full)],1);
    curve2_sort = sortrows([angle(curve2_full) abs(curve2_full)],1);
    
    % pad sorted curves with endpoints to wrap around
    curve1_sort_padded = [curve1_sort(end,:);curve1_sort;curve1_sort(1,:)];
    curve2_sort_padded = [curve2_sort(end,:);curve2_sort;curve2_sort(2,:)];
    curve1_sort_padded(1,1) = curve1_sort_padded(1,1)-2*pi;
    curve2_sort_padded(1,1) = curve2_sort_padded(1,1)-2*pi;
    curve1_sort_padded(end,1) = curve1_sort_padded(end,1)+2*pi;
    curve2_sort_padded(end,1) = curve2_sort_padded(end,1)+2*pi;
    
    % find distance at many rays
    rays = linspace(-pi,pi,10000)';
    rays = rays(2:end);
    curve1_mags = interp1(curve1_sort_padded(:,1),curve1_sort_padded(:,2),rays);
    curve2_mags = interp1(curve2_sort_padded(:,1),curve2_sort_padded(:,2),rays);
    dist = mean((curve1_mags-curve2_mags).^2);
end

function [rays,mags] = get_full_curve(curve)
    % pad curves with wrap-around term for full curve
    curve_wrap = [curve(end);curve];
    
    % find full curves and remove duplicate point
    interpolater = linspace(1,length(curve_wrap),length(curve_wrap)*500+1)';
    interpolater = interpolater(2:end);
    curve_full = interp1(curve_wrap,interpolater);
    
    % get sorted curves in angles and magnitudes
    curve_sort = sortrows([angle(curve_full) abs(curve_full)],1);
    
    % pad sorted curves with endpoints to wrap around
    curve_sort_padded = [curve_sort(end,:);curve_sort;curve_sort(1,:)];
    curve_sort_padded(1,1) = curve_sort_padded(1,1)-2*pi;
    curve_sort_padded(end,1) = curve_sort_padded(end,1)+2*pi;
    
    % find distance at many rays
    rays = linspace(-pi,pi,9)';
%     rays = rays(2:end);
    mags = interp1(curve_sort_padded(:,1),curve_sort_padded(:,2),rays);
end