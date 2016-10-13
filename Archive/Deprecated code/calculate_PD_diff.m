function dPD_data = calculate_PD_diff(bdf_cell_1,bdf_cell_2,leave_sorted,plotPDs)
% CALCULATE_PD_DIFF Finds the difference in preferred directions of units
% between two data sets.
%   CALCULATE_PD_DIFF(bdf_cell_1,bdf_cell_2,leave_sorted,plotPDs) uses
%   bootstrapping to calculate the difference in preferred direction for
%   units between bdf_cell_1 and bdf_cell_2. leave_sorted=1 for single
%   unit analysis, =0 for multi-unit analysis. BDFs should have
%   been merged and sorted together. plotPDs is an boolean indicating
%   whether the PDs should be plotted or not (1=yes, 0=no).
%
%   dPD_data output is a <number of units> x 7 array of dPD data. Columns
%   are:
%   [Channel #, Unit #, dPD (deg), Mod. depth, Low CI, High CI, CI range]

% Author: Raeed Chowdhury
% Date Revised: 07/08/2014

if(~leave_sorted)
    for i=1:length(bdf_cell_1)
        bdf_cell_new_1{i} = remove_sorting(bdf_cell_1{i});
    end
    for i=1:length(bdf_cell_2)
        bdf_cell_new_2{i} = remove_sorting(bdf_cell_2{i});
    end
else
    bdf_cell_new_1 = bdf_cell_1;
    bdf_cell_new_2 = bdf_cell_2;
end

% find valid units
ul_1 = unit_list(bdf_cell_new_1{1},~leave_sorted);
ul_2 = unit_list(bdf_cell_new_2{1},~leave_sorted);
ul = intersect(ul_1,ul_2,'rows');

ts = 50;
[pos,vel,acc] = sample_bdf_kin(bdf_cell_new_1,ts);
glm_input = [pos vel sqrt(vel(:,1).^2+vel(:,2).^2)];

% loop over units
dPD_data = zeros(length(ul),7);
tic;
for i = 1:length(ul)
    et = toc;
    fprintf(1, 'ET: %f (%d of %d)\n', et, i, length(ul));
    
    % Set up model inputs
    s_1 = bin_spikes(bdf_cell_new_1,ts,ul(i,1),ul(i,2));

    % find indices to keep
    keep_idx = 1:length(s_1);
%     glm_input = glm_input(keep_idx,:);
%     s = s(keep_idx);

    % GLM fit both workspaces
    b_bootstrap = zeros(6,1000);
    PD_boot = zeros(1,1000);
    for repct = 1:1000
        idx_rand = uint32(1+(length(keep_idx)-1)*rand(length(keep_idx),1));
        b_bootstrap(:,repct) = glmfit(glm_input(idx_rand,:),s_1(idx_rand),'poisson');
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
    
    dPD_data(i,:) = [double(ul(i,1)) double(ul(i,2)) PD_simple(i) moddepth(i) CI_boot(1,i) CI_boot(2,i) CI_range(i)];
end

if(plotPDs)
    figure
    for i=1:length(ul)
        subplot(4,8,i)
        polar(pi,1,'.')
        hold on
        polar([PD_simple(i) PD_simple(i)]*pi/180,[0 moddepth(i)/mean(moddepth)],'b')

        [x1,y1]=pol2cart(PD_simple(i)*pi/180,moddepth(i)/mean(moddepth)); % needed to fill up the space between the two CI
        [x2,y2]=pol2cart(CI_boot(1,i)*pi/180,moddepth(i)/mean(moddepth));
        [x3,y3]=pol2cart(CI_boot(2,i)*pi/180,moddepth(i)/mean(moddepth));
        x_fill = [0 x2 x1 x3];
        y_fill = [0 y2 y1 y3];
        patch(x_fill,y_fill,'b','facealpha',0.3);

        title(['Chan ' num2str(ul(i,1)) ', Unit ' num2str(ul(i,2))])
    end
end
